// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

//--------------------------------------------------------------------
//
// Interface to Trilinos linear solver library.
//
//--------------------------------------------------------------------

/*!
  \file    trilinos_linear_solver.cpp
  \brief   wrap Trilinos solver functions
*/

#include <cstdlib>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include "trilinos_impl.h"
#include "Profiler.h"
#include "ComMod.h"
#define NOOUTPUT

// The per-equation linear system state that used to live here -- dof,
// ghostAndLocalNodes, localNodes, coupledBC and the index vectors -- now lives
// in struct Trilinos (trilinos_impl.h). It was only safe as file-scope state
// because trilinos_lhs_create rebuilt it before every solve; now that the
// structure is reused across solves, a multi-physics run would otherwise solve
// one equation using another equation's dof and column indices.

int timecount = 0;

// ----------------------------------------------------------------------------
/**
 * Define the matrix vector multiplication operation to do at each iteration of
 * an iterative linear solver. This function is called by the AztecOO
 * (K + v1*v1' + v2*v2' + ...) *x = K*x + v1*(v1'*x) + v2*(v2'*x), 
 * where [v1, v2, ...] = bdryVec_list
 *
 * For coupled Neumann boundary terms (v1*v1', v2*v2'), we use efficient an vectorized 
 * operation rather than explicitly forming the rank 1 outer product matrix and 
 * adding it to the global stiffness matrix K.
 *
 * \param x vector to be applied on the operator
 * \param y result of sprase matrix vector multiplication
 */
void TrilinosMatVec::apply(const Tpetra_MultiVector& x, Tpetra_MultiVector& y,
           Teuchos::ETransp mode, Scalar_d alpha, Scalar_d beta) const
{
  // Store initial matrix vector product result in y (y = K*x if alpha = 1.0, beta = 0.0) 
  trilinos_->K->apply(x, y, mode, alpha, beta); //Y = beta*Y + alpha*K*X 

  // Now add on the coupled Neumann boundary contribution y += v1*(v1'*x) + v2*(v2'*x) + ...
  if (trilinos_->coupledBC)
  {
    // Cap terms contribute to the flux scalar term only, not to pressure update rows.
    // So we compute (v_face^T x + v_cap^T x) and apply it with v_face.
    for (size_t i = 0; i < trilinos_->bdryVec_list.size(); ++i)
    {
      Scalar_d dot_face = 0.0;
      Scalar_d dot_cap = 0.0;
      Teuchos::Array<Scalar_d> dots(1);

      auto bdryVec = trilinos_->bdryVec_list[i];
      x.dot(*bdryVec, dots());
      dot_face = dots[0];

      if (i < trilinos_->bdryCapVec_list.size()) {
        auto bdryCapVec = trilinos_->bdryCapVec_list[i];
        x.dot(*bdryCapVec, dots());
        dot_cap = dots[0];
      }

      y.update(dot_face + dot_cap, *bdryVec, 1.0);
    }
  }
}

// ----------------------------------------------------------------------------
/**
 * make the graph global and optimize storage on it
 * this function creates the LHS structure topology and RHS vector based on the
 * block map
 * \param numGlobalNodes total/global number of nodes each node can have dof
 *                       for the spatial dim
 * \param numLocalNodes  number of nodes owned by this proc in block coordiantes
 * \param numGhostAndLocalNodes number of nodes owned and shared by this
 *                              processor includes ghost nodes
 * \param nnz            number of nonzeros in LHS matrix of calling proc
 * \param ltgSorted      integer pointer of size numLocalNodes returns global
 *                       node number of local node a to give blockrow
 * \param ltgUnsorted    unsorted/not-reordered version
 * \param rowPtr         CSR row ptr of size numLocalNodes + 1 to block rows
 * \param colInd         CSR column indices ptr (size nnz points) to block rows
 * \param Dof            size of each block element to give dim of each block
 * 
 * \param numCoupledNeumannBC number of coupled Neumann BCs

 */

/// @brief Order-sensitive hash of an index vector, for the structure signature.
static std::size_t hash_index_vector(const Vector<int>& v)
{
  std::size_t h = 1469598103934665603ull;          // FNV-1a offset basis
  const std::size_t prime = 1099511628211ull;
  const int n = v.size();
  h = (h ^ static_cast<std::size_t>(n)) * prime;
  for (int i = 0; i < n; i++) {
    h = (h ^ static_cast<std::size_t>(static_cast<unsigned int>(v(i)))) * prime;
  }
  return h;
}

void trilinos_lhs_create(const Teuchos::RCP<Trilinos> &trilinos_, const int numGlobalNodes, const int numLocalNodes,
        const int numGhostAndLocalNodes, const int nnz, const Vector<int>& ltgSorted,
        const Vector<int>& ltgUnsorted, const Vector<int>& rowPtr, const Vector<int>& colInd,
        const int Dof, const int cpp_index, const int proc_id, const int numCoupledNeumannBC,
        const Array<double>& nodeCoords)
{
  // Everything built below is a function of these inputs alone, and they are
  // fixed for a given equation and mesh -- the sparsity pattern comes from
  // com_mod.rowPtr/colPtr, written only by lhsa() on the initialize path. This
  // used to be rebuilt on every Newton iteration. Rebuild only when something
  // actually changed: a different equation, or a remesh (which re-runs
  // initialize and so changes the signature).
  LhsSignature sig;
  sig.gtnNo = numGlobalNodes;
  sig.mynNo = numLocalNodes;
  sig.tnNo = numGhostAndLocalNodes;
  sig.nnz = nnz;
  sig.dof = Dof;
  sig.numCoupledNeumannBC = numCoupledNeumannBC;
  sig.ltg_hash = hash_index_vector(ltgSorted) ^ (hash_index_vector(ltgUnsorted) << 1);
  sig.rowPtr_hash = hash_index_vector(rowPtr);
  sig.colPtr_hash = hash_index_vector(colInd);

  if (trilinos_->signature.valid() && trilinos_->signature == sig) {
    // Structure unchanged, so the Map, ghostMap, graph, importer and vectors are
    // all still valid and the AMG hierarchy is still applicable.
    //
    // Only the matrix is recreated. It is released at the end of each solve
    // because reusing it changes the order in which values are summed into the
    // CSR (pre- versus post-fillComplete paths) and perturbs results at ~7e-10.
    // Rebuilding it from the already fill-complete graph is the cheap part.
    if (trilinos_->K.is_null()) {
      trilinos_->K = Teuchos::rcp(new Tpetra_CrsMatrix(trilinos_->K_graph));
    }
    return;
  }
  trilinos_->signature = sig;

  // Anything built from the old structure is invalid, including the AMG
  // hierarchy, whose grid transfers are tied to the graph.
  trilinos_->MueluPrec = Teuchos::null;
  trilinos_->ifpackPrec = Teuchos::null;

  #ifdef debug_trilinos_lhs_create
  std::string msg_prefix;
  msg_prefix = std::string("[trilinos_lhs_create:") + std::to_string(proc_id) + "] ";
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== trilinos_lhs_create ==========" << std::endl;
  std::cout << msg_prefix << "Dof: " << Dof << std::endl;
  std::cout << msg_prefix << "cpp_index: " << cpp_index << std::endl;
  #endif

  int indexBase = 1; //0 based indexing for C/C++ / 1 for Fortran
  if (cpp_index == 1) {
    indexBase = 0; 
  }

  trilinos_->dof = Dof; //constant size trilinos_->dof blocks
  trilinos_->ghostAndLocalNodes = numGhostAndLocalNodes;
  trilinos_->localNodes = numLocalNodes;

  trilinos_->comm = Tpetra::getDefaultComm();

  #ifdef debug_trilinos_lhs_create
  std::cout <<  msg_prefix << "indexBase: " << indexBase << std::endl;
  std::cout << msg_prefix << "trilinos_->dof: " << trilinos_->dof << std::endl;
  std::cout << msg_prefix << "trilinos_->ghostAndLocalNodes: " << trilinos_->ghostAndLocalNodes << std::endl;
  std::cout << msg_prefix << "trilinos_->localNodes: " << trilinos_->localNodes << std::endl;
  std::cout << msg_prefix << "trilinos_->localToGlobalSorted.size(): " << trilinos_->localToGlobalSorted.size() << std::endl;
  #endif

  // The variables need to be reset at every non-linear iterations
  // Tpetra is more strict then Epetra in graph allocation and once 
  // created does not allow any changes unless destroyed and reallocated
  trilinos_->localToGlobalSorted.clear(); 
  trilinos_->localToGlobalUnsorted.clear();
  trilinos_->globalColInd.clear();

  trilinos_->globalDofGIDs.clear();
  trilinos_->globalGhostDofGIDs.clear();
  

  // allocate memory for vectors
  trilinos_->localToGlobalSorted.reserve(numGhostAndLocalNodes);
  trilinos_->localToGlobalUnsorted.reserve(numGhostAndLocalNodes);
  trilinos_->globalColInd.reserve(nnz);

  trilinos_->globalDofGIDs.reserve(numGhostAndLocalNodes * trilinos_->dof);
  trilinos_->globalGhostDofGIDs.reserve(numGhostAndLocalNodes * trilinos_->dof); 

  trilinos_->nnzPerRow.clear (); // for fsils assembly
  trilinos_->nnzPerRow.reserve(numGhostAndLocalNodes);

  // Define localtoglobal to be used for unqiue partition map
  //only take ltgSorted(1:numLocalNodes) since those are owned by the processor
  //
  for (unsigned i = 0; i < numLocalNodes; ++i)
  {
    // any nodes following are ghost nodes so that subset
    for (int d = 0; d < trilinos_->dof; ++d)
    {
      trilinos_->globalDofGIDs.emplace_back(ltgSorted[i] * trilinos_->dof + d);
    }
  }

  for (unsigned i = 0; i < numGhostAndLocalNodes; ++i)
  {
    trilinos_->localToGlobalSorted.emplace_back(ltgSorted[i]);
    trilinos_->localToGlobalUnsorted.emplace_back(ltgUnsorted[i]);
    for (int d = 0; d < trilinos_->dof; ++d)
    {
      trilinos_->globalGhostDofGIDs.emplace_back(ltgSorted[i] * trilinos_->dof + d);
    }
    trilinos_->nnzPerRow.emplace_back(rowPtr[i+1] - rowPtr[i]);
  }

  /*
    Creating a Map for the local nodes owned by the processor
  */
  // Create the Tpetra Map (DoF-wise map)
  trilinos_->Map = Teuchos::rcp(new Tpetra_Map(
    Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
    Teuchos::arrayView(trilinos_->globalDofGIDs.data(), trilinos_->globalDofGIDs.size()), indexBase, trilinos_->comm));
  // Node-wise counterpart of Map: the same owned nodes without the dof blocking,
  // which is the layout MueLu expects for coordinates.
  const int nsd = nodeCoords.nrows();
  trilinos_->nodeMap = Teuchos::rcp(new Tpetra_Map(
    Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
    Teuchos::arrayView(ltgSorted.data(), numLocalNodes), indexBase, trilinos_->comm));

  {
    trilinos_->coords = Teuchos::rcp(new Tpetra_MultiVector(trilinos_->nodeMap, nsd));
    auto view = trilinos_->coords->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (int i = 0; i < numLocalNodes; i++) {
      for (int j = 0; j < nsd; j++) {
        view(i, j) = nodeCoords(j, i);
      }
    }
  }

  // Rigid body modes on the dof map. Translations are unit vectors in each
  // component; the rotations are the linearisation of a rotation about each
  // axis, so for a node at (x,y,z) they are (-y,x,0), (0,-z,y) and (z,0,-x).
  // Coordinates are shifted by their mean first, which keeps the rotational
  // modes well scaled against the translations on meshes sited far from the
  // origin. Only offered to MueLu when SVMP_MUELU_NULLSPACE=rbm.
  if (Dof >= nsd) {
    const int n_modes = (nsd == 3) ? 6 : 3;
    trilinos_->nullspace = Teuchos::rcp(new Tpetra_MultiVector(trilinos_->Map, n_modes));
    trilinos_->nullspace->putScalar(0.0);

    std::vector<double> centre(nsd, 0.0);
    for (int j = 0; j < nsd; j++) {
      double local_sum = 0.0;
      for (int i = 0; i < numLocalNodes; i++) { local_sum += nodeCoords(j, i); }
      double global_sum = 0.0;
      long local_n = numLocalNodes, global_n = 0;
      Teuchos::reduceAll(*trilinos_->comm, Teuchos::REDUCE_SUM, 1, &local_sum, &global_sum);
      Teuchos::reduceAll(*trilinos_->comm, Teuchos::REDUCE_SUM, 1, &local_n, &global_n);
      centre[j] = (global_n > 0) ? global_sum / static_cast<double>(global_n) : 0.0;
    }

    auto ns = trilinos_->nullspace->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (int i = 0; i < numLocalNodes; i++) {
      const double x = nodeCoords(0, i) - centre[0];
      const double y = nodeCoords(1, i) - centre[1];
      const double z = (nsd == 3) ? nodeCoords(2, i) - centre[2] : 0.0;
      for (int d = 0; d < nsd; d++) { ns(i * Dof + d, d) = 1.0; }
      if (nsd == 3) {
        ns(i * Dof + 0, 3) = -y;  ns(i * Dof + 1, 3) =  x;
        ns(i * Dof + 1, 4) = -z;  ns(i * Dof + 2, 4) =  y;
        ns(i * Dof + 0, 5) =  z;  ns(i * Dof + 2, 5) = -x;
      } else {
        ns(i * Dof + 0, 2) = -y;  ns(i * Dof + 1, 2) =  x;
      }
    }
  }

  /*
    Creating a Map for the local nodes owned and shared by the processor
  */
  // Create the ghost map — include owned + ghost GIDs
  trilinos_->ghostMap = Teuchos::rcp(new Tpetra_Map(
    Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
    Teuchos::arrayView(trilinos_->globalGhostDofGIDs.data(), trilinos_->globalGhostDofGIDs.size()), indexBase, trilinos_->comm));

  /*
    Graph construction  
  */
  // Calculate nnzPerDofRow to pass into graph constructor
  std::unordered_map<GO, size_t> gidToUnsortedIndex;
  gidToUnsortedIndex.reserve(ltgUnsorted.size());
  for (size_t i = 0; i < ltgUnsorted.size(); ++i)
    gidToUnsortedIndex[ltgUnsorted[i]] = i;

  const LO numLocalDofs = trilinos_->Map->getLocalNumElements();
  trilinos_->nnzPerDofRow.clear();
  trilinos_->nnzPerDofRow.reserve(numLocalDofs);

  for (LO lid = 0; lid < numLocalDofs; ++lid)
  {
    GO dofGid = trilinos_->Map->getGlobalElement(lid);
    GO nodeGid = dofGid / trilinos_->dof; 
    int d = dofGid % trilinos_->dof;

    auto it = gidToUnsortedIndex.find(nodeGid);
    if (it == gidToUnsortedIndex.end())
    {
      std::cerr << "[ERROR] nodeGid " << nodeGid << " not found in ltgUnsorted\n";
      continue;
    }

    size_t unsortedIdx = it->second;
    int numNodeCols = rowPtr[unsortedIdx + 1] - rowPtr[unsortedIdx];
    trilinos_->nnzPerDofRow.emplace_back(numNodeCols * trilinos_->dof);
  }

  // Construct graph based on nnz per row
  //
  trilinos_->K_graph = Teuchos::rcp(new Tpetra_CrsGraph(trilinos_->Map, trilinos_->nnzPerDofRow));

  unsigned nnzCount = 0; //cumulate count of block nnz per rows

  // Use unsortedltg to map local col ind to global indices
  //
  if (trilinos_->globalColInd.size() != nnz) { // only if nnz changed
    trilinos_->globalColInd.clear(); // destroy
    for (unsigned i = 0; i < nnz; ++i) {
      // Convert to global indexing subtract 1 for fortran vs C array indexing
      trilinos_->globalColInd.emplace_back(ltgUnsorted[colInd[i] - indexBase]);
    }
  } 
  
  //loop over block rows owned by current proc using localToGlobal index pointer
  for (size_t i = 0; i < numGhostAndLocalNodes; ++i) {
      GO nodeGID = ltgUnsorted[i];
      int numEntries = rowPtr[i+1] - rowPtr[i];

      for (int d = 0; d < trilinos_->dof; ++d) {
          GO rowGID = nodeGID * trilinos_->dof + d;
          std::vector<GO> rowCols(numEntries*trilinos_->dof);

          for (int j = 0; j < numEntries; ++j) {
              GO colNode = trilinos_->globalColInd[nnzCount + j];
              for (int col_d = 0; col_d < trilinos_->dof; ++col_d)
                  rowCols[j*trilinos_->dof + col_d] = colNode * trilinos_->dof + col_d;
          }

          trilinos_->K_graph->insertGlobalIndices(rowGID, rowCols);
      }
      nnzCount += numEntries;
  }

  // by end of iterations nnzCount should equal nnz-otherwise there is an error
  // Dofs are not counted in nnzCount, so we do not multiply by dof
  if (nnzCount != nnz)
  {
    std::cout << "Number of Entries Entered in Graph does not equal nnz " << std::endl;
    exit(1);
  }

  //check if trilinos methods are successful
  trilinos_->K_graph->fillComplete();
  if (!trilinos_->K_graph->isFillComplete())
  {
    std::cout << "ERROR: fillComplete() did not succeed on Graph" << std::endl;
    exit(1);
  }

  // --- Create block finite element matrix from graph with fillcomplete ------
  // construct matrix from filled graph
  trilinos_->K = Teuchos::rcp(new Tpetra_CrsMatrix(trilinos_->K_graph));

  //construct RHS force vector F topology
  trilinos_->F = Teuchos::rcp(new Tpetra_MultiVector(trilinos_->Map, 1));

  //construct RHS force vector ghostF topology
  trilinos_->ghostF = Teuchos::rcp(new Tpetra_MultiVector(trilinos_->ghostMap, 1));

  // Construct a boundary vector for each coupled Neumann boundary condition
  trilinos_->bdryVec_list.clear();
  trilinos_->bdryCapVec_list.clear();
  for (int i = 0; i < numCoupledNeumannBC; ++i)
  {
    trilinos_->bdryVec_list.push_back(Teuchos::rcp(new Tpetra_MultiVector(trilinos_->Map, 1)));
    trilinos_->bdryCapVec_list.push_back(Teuchos::rcp(new Tpetra_MultiVector(trilinos_->Map, 1)));
  }

  // Initialize solution vector which is unique and does not include the ghost
  // indices using the unique map
  trilinos_->X = Teuchos::rcp(new Tpetra_Vector(trilinos_->Map));

  //initialize vector which will import the ghost nodes using the ghost map
  trilinos_->ghostX = Teuchos::rcp(new Tpetra_Vector(trilinos_->ghostMap));

  //Create importer of the two maps
  trilinos_->Importer = Teuchos::rcp(new Tpetra_Import(trilinos_->Map, trilinos_->ghostMap));

} // trilinos_lhs_create_

// ----------------------------------------------------------------------------
/**
 * This function assembles the dense element block stiffness matrix and element
 * block force vector into the global block stiffness K and global block force
 * vector F.
 * Happens on a single element and will be called repeatedly as we loop over
 * all the elements.
 *
 * \param numNodesPerElement number of nodes (d) for element e
 * \param eqN                converts local element node number to local proc
 *                           equation number-size is numNodesPerElement
 * \param lK                 element local stiffness of size
 *                           (dof*dof, numNodesPerElement, numNodesPerElement)
 * \param lR                 element local force vector of size
 *                           (dof, numNodesPerElement)
 */
void trilinos_doassem_(const Teuchos::RCP<Trilinos> &trilinos_, int &numNodesPerElement, const int *eqN, const double *lK, double *lR)
{
  #ifdef debug_trilinos_doassem
  std::cout << "[trilinos_doassem_] ========== trilinos_doassem_ ===========" << std::endl;
  std::cout << "[trilinos_doassem_] trilinos_->dof: " << trilinos_->dof << std::endl;
  std::cout << "[trilinos_doassem_] numNodesPerElement: " << numNodesPerElement << std::endl;
  #endif

  //dof values per global ID in the force vector
  int numValuesPerID = trilinos_->dof;

  //converts eqN in local proc values to global values
  std::vector<int> localToGlobal(numNodesPerElement);

  for (int i = 0; i < numNodesPerElement; ++i)
    localToGlobal[i] = trilinos_->localToGlobalUnsorted[eqN[i]];

  //loop over local nodes on the element
  for (int a = 0; a < numNodesPerElement; ++a)
  {
    // Sum into contributions from element node-global assemble 
    // we assemble owned and ghost nodes (to be assembled across 
    // processors later)
    GO globalRow = localToGlobal[a] * trilinos_->dof; // first trilinos_->dof of node a
    
    // Sum into each DoF for vector F
    for (int d = 0; d < trilinos_->dof; ++d)
    {
      trilinos_->ghostF->sumIntoGlobalValue(globalRow + d, 0, lR[a * trilinos_->dof + d]);
    }
    
    // For the matrix K update:
    // For each node b connected to node a, insert the block values
    for (int b = 0; b < numNodesPerElement; ++b)
    {
      // The global column base for node b (all its dofs)
      GO colBase = localToGlobal[b] * trilinos_->dof;

      // We need to build row indices and values for each DoF row of node a
      for (int i = 0; i < trilinos_->dof; ++i)
      {
        GO globalRowK = globalRow + i;

        std::vector<GO> cols(trilinos_->dof);
        std::vector<Scalar_d> vals(trilinos_->dof);

        for (int j = 0; j < trilinos_->dof; ++j)
        {
          cols[j] = colBase + j;
          // 'a' and 'b' have to be inverted since lK.data() order them in a transpose way
          vals[j] = lK[b * trilinos_->dof * trilinos_->dof * numNodesPerElement + a * trilinos_->dof * trilinos_->dof + i * trilinos_->dof + j];
        }
        // sumIntoGlobalValues performs local to processor assembly 
        trilinos_->K->sumIntoGlobalValues(globalRowK, trilinos_->dof, vals.data(), cols.data());
      }
    }
  }
} // trilinos_doassem_

// ----------------------------------------------------------------------------
/**
 * Tthis function creates the LHS structure topology and RHS vector based on
 * the block map using the global stiffness matrix and global force vector
 * assuming the assembly part has been done in svFSI
 *
 * \param Val         nonzero values of assembled global stiffness
 * \param RHS         global force vector
 * \param x           solution vector
 * \param dirW        gives information about the Dirichlet BC
 * \param resNorm     norm of solution
 * \param initNorm    initial norm of RHS
 * \param numIters    linear solver number of iterations
 * \param solverTime  time in linear solver
 * \param dB          log ratio for output
 * \param converged   neither or not converged in the max number of iterations
 * \param lsType      defines type of linear solver to use
 * \param relTol      default relative tolerance or as set in input file
 * \param maxIters    default max number of iterations for gmres per restart
 * \param kspace      specific for gmres dim of the stored Krylov space vectors
 * \param precondType defines type of preconditioner to use
 */

void trilinos_global_solve_(const Teuchos::RCP<Trilinos> &trilinos_, const double *Val, const double *RHS, double *x,
        const double *dirW, double &resNorm, double &initNorm, int &numIters,
        double &solverTime, double &dB, bool &converged, int &lsType,
        double &relTol, int &maxIters, int &kspace, int &precondType)
{
  // Copying the FSILS-assembled values into the Tpetra matrix. Independent of the
  // preconditioner -- no hierarchy exists yet -- but it is a sizeable share of
  // the solve, so it is timed rather than left unaccounted.
  {
  SVMP_PROFILE_PHASE("ls_fill");

  int nnzCount = 0; //cumulate count of block nnz per rows
  int count = 0;
  int numValuesPerID = trilinos_->dof; //trilinos_->dof values per id pointer to trilinos_->dof
  std::vector<Scalar_d> values(trilinos_->dof); // holds local matrix entries
  std::vector<GO> colGIDK(trilinos_->dof); // holds global column indices for K

  // loop over block rows owned by current proc using localToGlobal index pointer
  //
  for (int i = 0; i < trilinos_->ghostAndLocalNodes; ++i) {
    int numEntries = trilinos_->nnzPerRow[i]; //block per of entries per row
    const GO rowGID = trilinos_->localToGlobalUnsorted[i];
    const GO* colGIDs = &trilinos_->globalColInd[nnzCount];

    // Copy RHS values into ghostF 
    for (int j = 0; j < trilinos_->dof; ++j) {
        trilinos_->ghostF->replaceGlobalValue(rowGID * trilinos_->dof + j, 0, RHS[i * trilinos_->dof + j]);
    }
    // Tpetra assembly from FSILS assembly, one call per matrix row rather than one
    // per (block column, row): each call converts the row GID to a local index and
    // binary-searches that row's column list, so batching pays it once per row.
    //
    // even though Val contains already assembled values we use sumIntoGlobalValues
    // because Tpetra can only store rows owned by the process and NOT ghost rows
    // when assemblying with FillComplete across processors missing contributions happens.
    // In this way we ensure all contributions are accounted for across processors, 
    // without duplicating contributions.
    const int dof = trilinos_->dof;
    values.resize(numEntries * dof);
    colGIDK.resize(numEntries * dof);
    for (int l = 0; l < dof; ++l) {
      for (int j = 0; j < numEntries; ++j) {
        for (int m = 0; m < dof; ++m) {
          colGIDK[j*dof + m] = colGIDs[j] * dof + m;
          values[j*dof + m] = Val[(count + j)*dof*dof + l*dof + m];
        }
      }
      trilinos_->K->sumIntoGlobalValues(rowGID * dof + l, numEntries * dof,
                                        values.data(), colGIDK.data());
    }
    count += numEntries;

    nnzCount += numEntries;
  }
  }  // ls_fill

  // Call solver code which assembles K and F for shared processors
  bool flagFassem = false;

  trilinos_solve_(trilinos_, x, dirW, resNorm, initNorm, numIters,
          solverTime, dB, converged, lsType,
          relTol, maxIters, kspace, precondType, flagFassem);

} // trilinos_global_solve_

// ----------------------------------------------------------------------------
/**
 * This function uses the established data structure to solve the block linear
 * system and passes the solution vector back to fortran with the ghost nodes
 * filled in
 *
 * \param x           solution vector
 * \param dirW        gives information about Dirichlet BC
 * \param resNorm     norm of the final residual
 * \param initNorm    norm of initial resiual x_init = 0 ||b||
 * \param numIters    number of iterations in the linear solver
 * \param solverTime  time in the linear solver
 * \param dB          log ratio for output
 * \param converged   can pass in solver and preconditioner type too
 * \param lsType      defines type of linear solver to use
 * \param relTol      default relative tolerance or as set in input file
 * \param maxIters    default max number of iterations for gmres per restart
 * \param kspace      specific for gmres dim of the stored Krylov space vectors
 * \param precondType defines type of preconditioner to use
 * \param isFassem    determines if F is already assembled at ghost nodes
 */
void trilinos_solve_(const Teuchos::RCP<Trilinos> &trilinos_, double *x, const double *dirW, double &resNorm,
        double &initNorm, int &numIters, double &solverTime, double &dB,
        bool &converged, int &lsType, double &relTol, int &maxIters,
        int &kspace, int &precondType, bool &isFassem)
{
  #define n_debug_trilinos_solve
  #ifdef debug_trilinos_solve
  std::cout << "[trilinos_solve] ========== trilinos_solve ==========" << std::endl;
  std::cout << "[trilinos_solve] resNorm: " << resNorm << std::endl;
  std::cout << "[trilinos_solve] initNorm: " << initNorm << std::endl;
  std::cout << "[trilinos_solve] lsType: " << lsType << std::endl;
  std::cout << "[trilinos_solve] precondType: " << precondType << std::endl;
  std::cout << "[trilinos_solve] isFassem: " << isFassem << std::endl;
  #endif
  bool flagFassem = isFassem;

  // Already filled from graph so does not need to call fillcomplete
  // routine will sum in contributions from elements on shared nodes amongst
  // processors
  //
  // Outlives the ls_setup scope: the solution is unscaled with it after the solve.
  Teuchos::RCP<Tpetra_Vector> diagonal = Teuchos::rcp(new Tpetra_Vector(trilinos_->Map));
  {
  SVMP_PROFILE_PHASE("ls_setup");
  // With a static graph the structure cannot change, so this is needed only
  // once -- and fillComplete() throws if called on an already-complete matrix.
  if (!trilinos_->K->isFillComplete()) {
    trilinos_->K->fillComplete();
  }

  // Importer runs in reverse here: doExport with an Import plan wants the plan's
  // target map to be the source object's (ghostMap) and its source map to be
  // this one's (Map), which is exactly how Importer was built. Constructing an
  // Export instead rebuilt the whole communication plan on every solve.
  if (flagFassem) {
    trilinos_->F->doExport(*trilinos_->ghostF, *trilinos_->Importer, Tpetra::ADD);
  } else { // RHS when using fsils assembly is already assembled and communicated
           // correctly among processors. REPLACE allows to create the correct
           // RHS vector of owned nodes only (no ghosts nodes)
    trilinos_->F->doExport(*trilinos_->ghostF, *trilinos_->Importer, Tpetra::REPLACE);
  }

  // Construct Jacobi scaling vector which uses dirW to take the Dirichlet BC
  // into account
  //
  constructJacobiScaling(trilinos_, dirW, *diagonal);

  // Compute norm of preconditioned multivector F
  Teuchos::Array<double> norms(1);
  trilinos_->F->norm2(norms());
  initNorm = norms[0];
  }  // ls_setup

  Teuchos::RCP<TrilinosMatVec> K_bdry = Teuchos::rcp(new TrilinosMatVec(trilinos_));

  // Define Belos linear problem if v is 0 does standard matvec product with K
  /*
    K_bdry -> Tpetra operator to include resistance boundary terms in LHS
    X -> Tpetra vector for the solution
    F -> Tpetra vector for the RHS
  */
  auto BelosProblem = Teuchos::rcp(new Belos_LinearProblem(K_bdry, trilinos_->X, trilinos_->F));

  {
    SVMP_PROFILE_PHASE("ls_precond");
    setPreconditioner(trilinos_, precondType, BelosProblem);
  }

  bool set = BelosProblem->setProblem();
  if (!set) {
    std::cerr << "ERROR: Belos LinearProblem setup failed!" << std::endl;
    exit(1);
  }
  
  Teuchos::ParameterList belosParams;
  int verbosity = Belos::Errors 
                + Belos::Warnings   
                + Belos::TimingDetails   
                + Belos::StatusTestDetails  
                + Belos::IterationDetails;         

  int maxRestarts = int(maxIters / kspace)+1;

  belosParams.set("Verbosity", verbosity);
  belosParams.set("Output Frequency", 1);
  belosParams.set("Output Style", 1);  
  belosParams.set("Convergence Tolerance", relTol);
  belosParams.set("Maximum Iterations", maxIters); 
  belosParams.set("Maximum Restarts", maxRestarts);
  
  std::string solverType;

  #ifdef NOOUTPUT
    belosParams.set("Verbosity", Belos::Errors);
  #endif

  // Solver is GMRES by default
  if (lsType == TRILINOS_GMRES_SOLVER) {
    solverType = "Block GMRES";
    belosParams.set("Num Blocks", kspace); 
    belosParams.set("Orthogonalization", "DGKS"); // DGKS orthogonalization is the most general
                                                  // and is recommended for most problems
  }
  else if (lsType == TRILINOS_BICGSTAB_SOLVER) {
    solverType = "BiCGStab";
  }
  else if (lsType == TRILINOS_CG_SOLVER) {
    solverType = "Pseudoblock CG";
  }
  else {
    throw std::runtime_error("Unknown solver type requested");
  }
  
  //checkStatus to calculate residual norm
  Belos_SolverFactory factory;
  auto solverManager = factory.create(solverType, Teuchos::rcpFromRef(belosParams));
  solverManager->setProblem(BelosProblem);

  // Run the solver (solve() handles restarts internally)
  converged = false;

  Teuchos::Time timer("Belos Solve Timer");
  timer.start();

  Belos::ReturnType result = solverManager->solve();
  
  timer.stop();

  solverTime = timer.totalElapsedTime();
  SVMP_PROFILE_ADD("ls_iterate", solverTime);

  if (result == Belos::Converged) {
    converged = true;
  }

  // Get number of iterations performed
  numIters = solverManager->getNumIters();

  // Get relative residual norm from the convergence test
  // Compute residual r = b - A*x manually since Belos does not provide it directly
  auto r = Teuchos::rcp(new Tpetra_MultiVector(trilinos_->F->getMap(), trilinos_->F->getNumVectors()));
  K_bdry->apply(*trilinos_->X, *r);
  r->update(1.0, *trilinos_->F, -1.0); // r = b - A*x

  // Compute residual norm
  Teuchos::Array<double> normR(1);
  r->norm2(normR());
  resNorm = normR[0];

  // Compute relative residual norm using precomputed initNorm
  double relRes = resNorm / initNorm;

  // Convert to decibel scale
  dB = 10.0 * log10(relRes);

  //Right scaling so need to multiply x by diagonal
  // Undo the right scaling and scatter the solution back to ghosted layout.
  SVMP_PROFILE_PHASE("ls_scatter");
  trilinos_->X->elementWiseMultiply(1.0, *trilinos_->X, *diagonal, 0.0);

  //Fill ghost X with x communicating ghost nodes amongst processors
  trilinos_->ghostX->doImport(*trilinos_->X, *trilinos_->Importer, Tpetra::INSERT);

  auto localView = trilinos_->ghostX->getLocalViewHost(Tpetra::Access::ReadOnly);
  size_t localLength = trilinos_->ghostX->getLocalLength();
  if (localLength == 0) {
    std::cout << "ERROR: Extracting copy of solution vector!" << std::endl;
    exit(1);
  }
  for (size_t i = 0; i < localLength; ++i) {
    x[i] = localView(i, 0);
  }

  // Zero out residual and solution vectors, lhs and, if used,
  // any preconditioner objects
  trilinos_->ghostF->putScalar(0.0);
  trilinos_->F->putScalar(0.0);
  if (trilinos_->coupledBC) {
    for (auto bdryVec : trilinos_->bdryVec_list)
      bdryVec->putScalar(0.0);
    for (auto bdryCapVec : trilinos_->bdryCapVec_list)
      bdryCapVec->putScalar(0.0);

  }
  trilinos_->X->putScalar(0.0);

  // The MueLu hierarchy is kept and refreshed on the next solve; see
  // setMueLuPreconditioner. Ifpack2 preconditioners are still dropped -- they
  // are cheap relative to an AMG setup and hold references to matrix values.
  if (trilinos_->ifpackPrec != Teuchos::null) trilinos_->ifpackPrec = Teuchos::null;

  // The matrix is released but the Map, ghostMap, graph, importer and vectors
  // are kept (see trilinos_lhs_create's signature check), so only the CrsMatrix
  // is reconstructed -- from an already fill-complete graph, which is cheap.
  //
  // Reusing the matrix itself was measured and does NOT preserve results: on a
  // freshly constructed matrix sumIntoGlobalValues accumulates into a
  // pre-fillComplete structure, while on a reused fill-complete matrix it takes
  // the local-index path. Same values, different summation order, ~7e-10
  // relative divergence -- enough to fail the 1e-10 regression tolerance.
  trilinos_->K = Teuchos::null;

} // trilinos_solve_

// ----------------------------------------------------------------------------
void setPreconditioner(const Teuchos::RCP<Trilinos> &trilinos_, int precondType, 
  Teuchos::RCP<Belos_LinearProblem>& BelosProblem)
{
  if (precondType == TRILINOS_DIAGONAL_PRECONDITIONER ||
      precondType == NO_PRECONDITIONER) {
    BelosProblem->setLeftPrec(Teuchos::null);
    return;
  }

  // Create Ifpack2 preconditioner
  Ifpack2::Factory factory;
  Teuchos::ParameterList precParams;    // Parameter list setup

  if (precondType == TRILINOS_BLOCK_JACOBI_PRECONDITIONER) {
    checkDiagonalIsZero(trilinos_);
    trilinos_->ifpackPrec = factory.create<Tpetra_CrsMatrix>("RELAXATION", trilinos_->K);
    precParams.set("relaxation: type", "Jacobi");
    precParams.set("relaxation: sweeps", 1);

  } else if (precondType == TRILINOS_ILU_PRECONDITIONER) {
    checkDiagonalIsZero(trilinos_);
    trilinos_->ifpackPrec = factory.create<Tpetra_CrsMatrix>("SCHWARZ", trilinos_->K);
    precParams.set("schwarz: inner preconditioner name", "ILUT");
    precParams.set("schwarz: combine mode", "Add");
    precParams.set("schwarz: overlap level", 1);
    precParams.set("fact: level-of-fill", 0); // Classic ILU(0)
     precParams.set("fact: relax value", 0.0);

  } else if (precondType == TRILINOS_ILUT_PRECONDITIONER) {
    checkDiagonalIsZero(trilinos_);
    trilinos_->ifpackPrec = factory.create<Tpetra_CrsMatrix>("SCHWARZ", trilinos_->K);
    precParams.set("schwarz: overlap level", 1);
    precParams.set("schwarz: inner preconditioner name", "ILUT");
    precParams.set("schwarz: combine mode", "Add");
    precParams.set("fact: ilut level-of-fill", 2.0);
    precParams.set("fact: drop tolerance", 1e-2);
    precParams.set("fact: relax value", 0.0);

  } else if (precondType == TRILINOS_RILUK0_PRECONDITIONER) {
    checkDiagonalIsZero(trilinos_);
    trilinos_->ifpackPrec = factory.create<Tpetra_CrsMatrix>("SCHWARZ", trilinos_->K);
    precParams.set("schwarz: inner preconditioner name", "RILUK"); 
    precParams.set("fact: level-of-fill", 0);
    precParams.set("fact: drop tolerance", 0.0);                   
    precParams.set("fact: relax value", 0.0);

  } else if (precondType == TRILINOS_RILUK1_PRECONDITIONER) {
    checkDiagonalIsZero(trilinos_);
    trilinos_->ifpackPrec = factory.create<Tpetra_CrsMatrix>("SCHWARZ", trilinos_->K);
    precParams.set("schwarz: inner preconditioner name", "RILUK"); 
    precParams.set("fact: level-of-fill", 1); 
    precParams.set("fact: drop tolerance", 1e-3);
    precParams.set("fact: relax value", 0.0);

  } else if (precondType == TRILINOS_ML_PRECONDITIONER) {
    checkDiagonalIsZero(trilinos_);
    setMueLuPreconditioner(trilinos_->MueluPrec, trilinos_->K, trilinos_->dof,
                           trilinos_->coords, trilinos_->nullspace);
    BelosProblem->setLeftPrec(trilinos_->MueluPrec);
    return;
  } else {
    throw std::runtime_error("[ERROR Trilinos] Unsupported preconditioner type.");
  }

  trilinos_->ifpackPrec->setParameters(precParams);
  trilinos_->ifpackPrec->initialize();
  trilinos_->ifpackPrec->compute();

  BelosProblem->setLeftPrec(trilinos_->ifpackPrec);

} // setPreconditioner

// ----------------------------------------------------------------------------
/*
 * Tune parameters for ML preconditioner using MueLu package
 * For a complete guide, refer to the MueLu documentation:
 * https://trilinos.github.io/pdfs/mueluguide.pdf
 */
/// How much of the AMG hierarchy to carry over between solves.
///
/// The matrix values change every Newton iteration but the sparsity pattern does
/// not, so the aggregation and grid transfers -- which depend only on the graph --
/// need not be recomputed. Accepted by MueLu: "none", "S", "tP", "RP", "emin",
/// "RAP", "full", in increasing order of what is kept.
///
///   none  rebuild everything (the original behaviour)
///   S     reuse smoothers only
///   tP    also reuse the tentative prolongator (skips aggregation)
///   RP    also reuse R and P (skips prolongator smoothing); coarse operators
///         are still recomputed from the new values
///   RAP   also reuse the coarse-grid pattern
///   full  reuse everything
///
/// Edit here to try a level. Higher reuse cuts setup time but can weaken the
/// preconditioner, which shows up as more Krylov iterations -- watch ls_iterate,
/// not just ls_precond. An XML option can follow once a level is settled on.
static constexpr const char* MUELU_REUSE_TYPE = "RAP";

void setMueLuPreconditioner(Teuchos::RCP<MueLu_Preconditioner> &MueLuPrec,
                            const Teuchos::RCP<Tpetra_CrsMatrix> &A,
                            const int dof,
                            const Teuchos::RCP<Tpetra_MultiVector> &coords,
                            const Teuchos::RCP<Tpetra_MultiVector> &nullspace)
{
  // 2c: an existing hierarchy is refreshed with the new matrix values rather
  // than rebuilt. MueluPrec is cleared by trilinos_lhs_create whenever the
  // structure signature changes, so a non-null value here means the pattern is
  // unchanged and reuse is valid.
  if (!MueLuPrec.is_null()) {
    auto muelu_op = Teuchos::rcp_dynamic_cast<
        MueLu::TpetraOperator<Scalar_d, LO, GO, Node>>(MueLuPrec);
    if (!muelu_op.is_null()) {
      MueLu::ReuseTpetraPreconditioner(A, *muelu_op);
      return;
    }
  }

  // MueLuPrec is now a Tpetra::Operator that can be plug into BelosProblem
  std::string optionsFile = "mueluOptions.xml";

  // The following built-in parameters proved to be generally good for big FSI problems
  Teuchos::ParameterList mueluParams;

  // MEASUREMENT: "high" prints the level-by-level row/nnz table and the
  // smoother chosen per level.
  const char* vb_env = std::getenv("SVMP_MUELU_VERBOSITY");
  mueluParams.set("verbosity", vb_env ? vb_env : "none");
  mueluParams.set("max levels", 6);   // number of multigrid levels
  mueluParams.set("cycle type", "V"); // V-cycle is standard
  mueluParams.set("fuse prolongation and update", true);

  // Problem type
  mueluParams.set("problem: type", "unknown"); // FSI is generally nonsymmetric
  mueluParams.set("reuse: type", MUELU_REUSE_TYPE);
  mueluParams.set("number of equations", dof);   //dof for this equation

  // Aggregation
  mueluParams.set("aggregation: type", "uncoupled");
  mueluParams.set("aggregation: min agg size", 2);
  mueluParams.set("aggregation: max agg size", 8);
  mueluParams.set("aggregation: ordering", "natural");
  mueluParams.set("aggregation: drop scheme", "classical");
  mueluParams.set("aggregation: strength-of-connection: measure", "smoothed aggregation");
  mueluParams.set("aggregation: number of random vectors", 5);
  mueluParams.set("aggregation: number of times to pre or post smooth", 3);

  // Smoothed Aggregation-specific parameters
  mueluParams.set("sa: damping factor", 1.0);
  mueluParams.set("sa: use filtered matrix", false);

  // Smoother
  mueluParams.set("smoother: type", "RELAXATION");
  mueluParams.set("smoother: pre or post", "both");
  mueluParams.set("smoother: overlap", 0);

  Teuchos::ParameterList &smootherParams = mueluParams.sublist("smoother: params");
  smootherParams.set("relaxation: type", "Gauss-Seidel");
  smootherParams.set("relaxation: sweeps", 2);

  // Coarse solver
  mueluParams.set("coarse: type", "KLU"); // robust direct solver for FSI
  mueluParams.set("coarse: max size", 2000);

  // Create MueLu preconditioner from matrix and parameter list, as Tpetra::Operator
  // MueLu reads these from the "user data" sublist and converts them to Xpetra
  // itself. Without coordinates it cannot do distance-based aggregation or
  // geometric repartitioning. The near-nullspace is opt-in: MueLu's default is
  // one constant vector per equation, i.e. the translations only, and swapping
  // in the six rigid body modes changes the coarse grid for every run.
  const char* ns_env = std::getenv("SVMP_MUELU_NULLSPACE");
  const bool use_rbm = (ns_env != nullptr && std::string(ns_env) == "rbm");
  if (!coords.is_null()) {
    mueluParams.sublist("user data").set("Coordinates", coords);
  }
  if (use_rbm && !nullspace.is_null()) {
    mueluParams.sublist("user data").set("Nullspace", nullspace);
  }

  std::ifstream ifs(optionsFile.c_str());
  if (ifs.good())
  {
    try
    {
        // Read the file into a list rather than handing MueLu the path. The path
        // overload builds its own list, so it would silently drop everything set
        // above -- the coordinates and near-nullspace, which the file cannot
        // express at all, and the two run-time parameters reinstated below.
        Teuchos::RCP<Teuchos::ParameterList> fileParams =
            Teuchos::getParametersFromXmlFile(optionsFile);

        // Deployment wiring, not AMG tuning: an options file should not have to
        // restate either, and losing them is silent. Wrong equations-per-node
        // costs convergence; losing reuse rebuilds the hierarchy every solve.
        // Set only when absent, so a file that states them deliberately wins.
        if (!fileParams->isParameter("number of equations")) {
          fileParams->set("number of equations", dof);
        }
        if (!fileParams->isParameter("reuse: type")) {
          fileParams->set("reuse: type", MUELU_REUSE_TYPE);
        }
        if (!coords.is_null()) {
          fileParams->sublist("user data").set("Coordinates", coords);
        }
        if (use_rbm && !nullspace.is_null()) {
          fileParams->sublist("user data").set("Nullspace", nullspace);
        }

        MueLuPrec = MueLu::CreateTpetraPreconditioner(
            Teuchos::rcp_static_cast<Tpetra_Operator>(A), *fileParams);
    }
    catch (std::exception &e)
    {
      std::cerr << "[MueLu Warning]: failed to create MueLu from file '" << optionsFile
                << "': " << e.what() << "\n  Falling back to built-in parameters.\n";
      MueLuPrec = MueLu::CreateTpetraPreconditioner(
          Teuchos::rcp_static_cast<Tpetra_Operator>(A), mueluParams);
    }
  }
  else
  {
    MueLuPrec = MueLu::CreateTpetraPreconditioner(
        Teuchos::rcp_static_cast<Tpetra_Operator>(A), mueluParams);
  }

}

// ----------------------------------------------------------------------------
/**
 * This routine is to be used with preconditioners such as BlockJacobi/ILU/ILUT
 * which require 1s on the diagonal
 */
void checkDiagonalIsZero(const Teuchos::RCP<Trilinos> &trilinos_)
{
  Teuchos::RCP<const Tpetra_Map> rowMap = trilinos_->K->getRowMap();
  Tpetra_Vector diagonal(rowMap);
  trilinos_->K->getLocalDiagCopy(diagonal);
  bool isZeroDiag = false;
  auto diagData = diagonal.getLocalViewHost(Tpetra::Access::ReadWrite);
  for (size_t i = 0; i < diagData.extent(0); ++i)
  {
    if (diagData(i, 0) == 0.0)
    {
      diagData(i, 0) = 1.0;
      isZeroDiag = true;
    }
  }
  Tpetra::replaceDiagonalCrsMatrix(*trilinos_->K, diagonal);

} // void checkDiagonalIsZero()

// ----------------------------------------------------------------------------
/**
 * To be called within solve-output diagonal from here
 *
 * \param dirW    pass in array with Dirichlet boundary face nodes marked
 * \paramdiagonal diagonal scaling vector need to output to multiply solution by
 */
void constructJacobiScaling(const Teuchos::RCP<Trilinos> &trilinos_, const double *dirW, Tpetra_Vector& diagonal)
{
  // Set Dirichlet weights
  for (int i = 0; i < trilinos_->localNodes; ++i) {
    for (int j = 0; j < trilinos_->dof; ++j) {
      // Map is built from globalDofGIDs, filled as ltgSorted[i]*dof + j over the
      // same i and j in the same order, so this entry's local index is i*dof + j.
      // Looking the GID up cost two hash probes to rediscover that.
      diagonal.replaceLocalValue(i * trilinos_->dof + j, dirW[i * trilinos_->dof + j]);
    }
  }

  // Extract and modify diagonal of K
  Tpetra_Vector Kdiag(trilinos_->K->getRowMap());
  trilinos_->K->getLocalDiagCopy(Kdiag);

  auto KdiagView = Kdiag.getLocalViewHost(Tpetra::Access::ReadWrite);
  for (size_t i = 0; i < KdiagView.extent(0); ++i) {
    if (KdiagView(i, 0) == 0.0)
      KdiagView(i, 0) = 1.0;
    KdiagView(i, 0) = 1.0 / std::sqrt(std::abs(KdiagView(i, 0)));
  }

  // diagonal = diagonal * Kdiag (element-wise)
  diagonal.elementWiseMultiply(1.0, diagonal, Kdiag, 0.0);

  // Apply scaling to K and F
  trilinos_->K->leftScale(diagonal);
  trilinos_->F->elementWiseMultiply(1.0, diagonal, *trilinos_->F, 0.0);
  trilinos_->K->rightScale(diagonal);

  // Scale boundary vectors if coupledBC is set
  if (trilinos_->coupledBC) {
    for (auto bdryVec : trilinos_->bdryVec_list) {
      bdryVec->elementWiseMultiply(1.0, diagonal, *bdryVec, 0.0);
    }
    for (auto bdryCapVec : trilinos_->bdryCapVec_list) {
      bdryCapVec->elementWiseMultiply(1.0, diagonal, *bdryCapVec, 0.0);
    }
  }
} // void constructJacobiScaling()

// ----------------------------------------------------------------------------
/**
 * \param  v            coupled boundary vector
 * \param  isCoupledBC  determines if coupled resistance BC is turned on
 */
void trilinos_bc_create_(const Teuchos::RCP<Trilinos> &trilinos_,
  const std::vector<Array<double>> &v_list, const std::vector<Array<double>> &vcap_list, bool &isCoupledBC)
{
  // store as global to determine which matvec multiply to use in solver
  trilinos_->coupledBC = isCoupledBC;

  if (isCoupledBC)
  {
    for (int i = 0; i < trilinos_->ghostAndLocalNodes; ++i)
    {
      auto globalRow = static_cast<GO>(trilinos_->localToGlobalSorted[i]);

      for (int k = 0; k < v_list.size(); ++k)
      {
        const double* v = v_list[k].data();
        auto bdryVec = trilinos_->bdryVec_list[k];
        const double* vcap = vcap_list[k].data();
        auto bdryCapVec = trilinos_->bdryCapVec_list[k];
        for (int j = 0; j < trilinos_->dof; ++j)
        {
          bdryVec->replaceGlobalValue(globalRow * trilinos_->dof + j, 0, v[i * trilinos_->dof + j]);
          bdryCapVec->replaceGlobalValue(globalRow * trilinos_->dof + j, 0, vcap[i * trilinos_->dof + j]);
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------
/**
 * for debugging purposes here are routines to print the matrix and RHS vector
 */
void printMatrixToFile(const Teuchos::RCP<Trilinos> &trilinos_)
{
  std::ofstream Kfile("K.txt");
  Kfile << std::scientific;
  Kfile.precision(17);

  Teuchos::RCP<const Tpetra_CrsMatrix> K = trilinos_->K;
  Teuchos::RCP<const Tpetra_Map> rowMap = K->getRowMap();

  for (int i = 0; i < trilinos_->ghostAndLocalNodes; ++i)
  {
    int nodeGlobalID = trilinos_->localToGlobalUnsorted[i];

    // Loop over each DoF in the block row
    for (int d = 0; d < trilinos_->dof; ++d)
    {
      int rowGID = nodeGlobalID * trilinos_->dof + d;

      // Get number of entries in the row
      size_t maxEntries = K->getNumEntriesInGlobalRow(rowGID);
      Tpetra_CrsMatrix::nonconst_global_inds_host_view_type indices("indices", maxEntries);
      Tpetra_CrsMatrix::nonconst_values_host_view_type values("values", maxEntries);

      size_t numEntries = 0;

      // Now get the actual data
      K->getGlobalRowCopy(rowGID, indices, values, numEntries);

      for (size_t j = 0; j < numEntries; ++j)
      {
        int colGID = indices(j);
        double val = values(j);
        Kfile << "Row " << rowGID << " Col " << colGID << " Val " << val << std::endl;
      }
    }
  }

  Kfile.close();
}

void printRHSToFile(const Teuchos::RCP<Trilinos> &trilinos_)
{
  std::ofstream Ffile("F.txt");
  Ffile.precision(17);

  // Get the local length of the Tpetra::MultiVector (assumes single vector)
  Teuchos::RCP<const Tpetra_MultiVector> F = trilinos_->F;
  std::size_t localLength = F->getLocalLength();

  // Create a view of the local data
  Teuchos::ArrayRCP<const double> F_local = F->getData(0); // 0 = first vector in MultiVector

  // Print each value to file
  for (std::size_t i = 0; i < localLength; ++i)
    Ffile << F_local[i] << std::endl;

  Ffile.close();
}

// ----------------------------------------------------------------------------
/**
 */
void printSolutionToFile(const Teuchos::RCP<Trilinos> &trilinos_)
{
  std::ofstream Xfile("X.txt");
  Xfile.precision(17);

  // Get the local length of the Tpetra::MultiVector (assumes single vector)
  Teuchos::RCP<const Tpetra_MultiVector> X = trilinos_->X;
  std::size_t localLength = X->getLocalLength();

  // Create a view of the local data
  Teuchos::ArrayRCP<const double> X_local = X->getData(0); // 0 = first vector in MultiVector

  // Print each value to file
  for (std::size_t i = 0; i < localLength; ++i)
    Xfile << X_local[i] << std::endl;

  Xfile.close();
}

/////////////////////////////////////////////////////////////////
//                  T r i l i n o s I m p l                    //
/////////////////////////////////////////////////////////////////

//--------------
// TrilinosImpl 
//--------------
// The TrilinosImpl private class hides Trilinos data structures
// and functions.
//
class TrilinosLinearAlgebra::TrilinosImpl {
  public:
    TrilinosImpl():trilinos_(Teuchos::rcp(new Trilinos())) {}
    void alloc(ComMod& com_mod, eqType& lEq);
    void assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN,
        const Array3<double>& lK, const Array<double>& lR);
    void initialize(ComMod& com_mod);
    void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res);
    void solve_assembled(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res);
    void init_dir_and_coup_neu(ComMod& com_mod, const Vector<int>& incL, const Vector<double>& res);
    void set_preconditioner(consts::PreconditionerType preconditioner);
    void finalize();

    consts::PreconditionerType preconditioner_;

    /// @brief Local to global mapping
    Vector<int> ltg_;

    /// @brief Factor for Dirichlet BCs
    Array<double> W_;

    /// @brief Residual
    Array<double> R_;
  
  private:
    Teuchos::RCP<Trilinos> trilinos_;
};

/// @brief Allocate Trilinos arrays.
void TrilinosLinearAlgebra::TrilinosImpl::alloc(ComMod& com_mod, eqType& lEq) 
{
  int dof = com_mod.dof;
  int tnNo = com_mod.tnNo;
  int gtnNo = com_mod.gtnNo;
  auto& lhs = com_mod.lhs;
  #ifdef n_debug_alloc
  std::cout << "[TrilinosImpl.alloc] dof: " << dof << std::endl;
  std::cout << "[TrilinosImpl.alloc] tnNo: " << tnNo << std::endl;
  std::cout << "[TrilinosImpl.alloc] gtnNo: " << gtnNo << std::endl;
  std::cout << "[TrilinosImpl.alloc] ltg_.size(): " << ltg_.size() << std::endl;
  #endif

  if (W_.size() != 0) {
    W_.clear();
    R_.clear();
  }

  W_.resize(dof,tnNo); 
  R_.resize(dof,tnNo);

  int cpp_index = 1;
  int task_id = com_mod.cm.idcm();

  // Nodal coordinates for the owned rows, in the same order as ltg_. That order
  // is the solver's node numbering permuted by lhs.map(), the permutation ltg_
  // itself is built with, so the two stay consistent.
  const int nsd = com_mod.nsd;
  Array<double> nodeCoords(nsd, lhs.mynNo);
  for (int a = 0; a < tnNo; a++) {
    const int k = com_mod.lhs.map(a);
    if (k < lhs.mynNo) {
      for (int j = 0; j < nsd; j++) {
        nodeCoords(j, k) = com_mod.x(j, a);
      }
    }
  }

  trilinos_lhs_create(trilinos_, gtnNo, lhs.mynNo, tnNo, lhs.nnz, ltg_, com_mod.ltg, com_mod.rowPtr, 
      com_mod.colPtr, dof, cpp_index, task_id, com_mod.lhs.nFaces, nodeCoords);

  // The matrix is now reused across solves rather than rebuilt, so it has to be
  // cleared here -- it used to arrive as a freshly constructed (zeroed) matrix
  // and the assembly paths only ever sum into it. Without this, values would
  // accumulate across Newton iterations.
  //
  // Done once per iteration, before any assembly: trilinos_doassem_ is called
  // per element and must not clear anything.
  if (!trilinos_->K.is_null()) {
    trilinos_->K->setAllToScalar(0.0);
  }
}

/// @brief Assemble local element arrays.
void TrilinosLinearAlgebra::TrilinosImpl::assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN,
        const Array3<double>& lK, const Array<double>& lR)
{
  trilinos_doassem_(trilinos_, const_cast<int&>(num_elem_nodes), eqN.data(), lK.data(), lR.data());
}

/// @brief Set data for Dirichlet and coupled Neumann boundary conditions.
void TrilinosLinearAlgebra::TrilinosImpl::init_dir_and_coup_neu(ComMod& com_mod, const Vector<int>& incL, const Vector<double>& res)
{
  using namespace consts;
  using namespace fsi_linear_solver;

  int dof = com_mod.dof;
  int gtnNo = com_mod.gtnNo;
  int tnNo = com_mod.tnNo;
  auto& lhs = com_mod.lhs;

  if (lhs.nFaces != 0) {
    for (auto& face : lhs.face) {
      face.incFlag = true;
    }

    for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
      if (incL(faIn) == 0)  {
        lhs.face[faIn].incFlag = false;
      }
    }

    for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
      auto& face = lhs.face[faIn];
      face.coupledFlag = false;
      if (!face.incFlag) {
        continue;
      }

      bool flag = (face.bGrp == BcType::BC_TYPE_Neu);
      if (flag && res(faIn) != 0.0) {
        face.res = res(faIn);
        face.coupledFlag = true;
      }
    }
  }

  W_ = 1.0;

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];
    if (!face.incFlag) {
      continue;
    }

    int faDof = std::min(face.dof,dof);

    if (face.bGrp == BcType::BC_TYPE_Dir) {
      for (int a = 0; a < face.nNo; a++) {
        int Ac = face.glob(a);
        for (int i = 0; i < faDof; i++) {
          W_(i,Ac) = W_(i,Ac) * face.val(i,a);
        }
      }
    }
  }

  std::vector<Array<double>> v_list;
  std::vector<Array<double>> vcap_list;
  bool isCoupledBC = false;

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    // Extract face
    auto& face = lhs.face[faIn]; 

    // Create a new array for each face and add it to the list
    v_list.push_back(Array<double>(dof,tnNo));
    auto& v = v_list.back();
    vcap_list.push_back(Array<double>(dof,tnNo));
    auto& vcap = vcap_list.back();
    
    if (face.coupledFlag) {
      isCoupledBC = true;
      int faDof = std::min(face.dof,dof);

      // Compute the coupled Neumann BC vector and store it in v (main face contribution)
      for (int a = 0; a < face.nNo; a++) {
        int Ac = face.glob(a);
        for (int i = 0; i < faDof; i++) {
          v(i,Ac) = v(i,Ac) + sqrt(fabs(res(faIn))) * face.val(i,a);
        }
      }
      
      // Add cap contribution (if cap exists for this face)
      if (face.cap_val.size() != 0 && face.cap_glob.size() != 0) {
        int cap_nNo = face.cap_val.ncols();
        for (int a = 0; a < cap_nNo; a++) {
          int Ac = face.cap_glob(a);
          if (Ac < 0) continue;  // cap node not on this rank
          for (int i = 0; i < faDof; i++) {
            vcap(i,Ac) = vcap(i,Ac) + sqrt(fabs(res(faIn))) * face.cap_val(i,a);
          }
        }
      }
    }
  }

  // Add the v vectors to global bdryVec_list
  trilinos_bc_create_(trilinos_, v_list, vcap_list, isCoupledBC);

}

/// @brief Initialze an array used for something.
void TrilinosLinearAlgebra::TrilinosImpl::initialize(ComMod& com_mod)
{
  #ifdef n_debug_initialize
  std::cout << "[TrilinosImpl] ---------- initialize ---------- " << std::endl;
  std::cout << "[TrilinosImpl.initialize] com_mod.tnNo: " << com_mod.tnNo << std::endl;
  #endif

  if (!Kokkos::is_initialized()) {
    Kokkos::initialize();
  }

  ltg_.resize(com_mod.tnNo);

  for (int a = 0; a < com_mod.tnNo; a++) {
    ltg_(com_mod.lhs.map(a)) = com_mod.ltg(a);
  }
}

void TrilinosLinearAlgebra::TrilinosImpl::finalize()
{
  // The structure and the AMG hierarchy are now kept between solves rather than
  // destroyed at the end of each one, so they are still alive here. They must be
  // released before Kokkos::finalize(), or Tpetra reports execution space
  // instances surviving past shutdown.
  if (!trilinos_.is_null()) {
    trilinos_->MueluPrec = Teuchos::null;
    trilinos_->ifpackPrec = Teuchos::null;
    trilinos_->K = Teuchos::null;
    trilinos_->K_graph = Teuchos::null;
    trilinos_->F = Teuchos::null;
    trilinos_->ghostF = Teuchos::null;
    trilinos_->X = Teuchos::null;
    trilinos_->ghostX = Teuchos::null;
    trilinos_->Importer = Teuchos::null;
    trilinos_->bdryVec_list.clear();
    trilinos_->bdryCapVec_list.clear();
    trilinos_->Map = Teuchos::null;
    trilinos_->ghostMap = Teuchos::null;
    trilinos_->comm = Teuchos::null;
    trilinos_->signature = LhsSignature();
  }

  if (Kokkos::is_initialized())
  {
    Kokkos::finalize();
  }
}

/// @brief Set the preconditioner.
void TrilinosLinearAlgebra::TrilinosImpl::set_preconditioner(consts::PreconditionerType prec_type)
{
  preconditioner_ = prec_type;
}

/// @brief Solve a system of linear equations assembled by fsils.
void TrilinosLinearAlgebra::TrilinosImpl::solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, 
    const Vector<double>& res)
{
  init_dir_and_coup_neu(com_mod, incL, res);

  auto& Val = com_mod.Val;
  auto& R = com_mod.R;
  int solver_type = static_cast<int>(lEq.ls.LS_type);
  int prec_type = static_cast<int>(preconditioner_);
  #define n_debug_solve
  #ifdef debug_solve
  std::cout << "[TrilinosImpl::solve] ---------- solve ---------- " << std::endl;
  auto prec_name = consts::preconditioner_type_to_name.at(preconditioner_); 
  std::cout << "[TrilinosImpl::solve] solver_type: " << solver_type << std::endl;
  std::cout << "[TrilinosImpl::solve] prec_type: " << prec_name << std::endl;
  std::cout << "[TrilinosImpl::solve] Val.size(): " << Val.size() << std::endl;
  std::cout << "[TrilinosImpl::solve] R_.size(): " << R_.size() << std::endl;
  std::cout << "[TrilinosImpl::solve] W_.size(): " << W_.size() << std::endl;
  std::cout << "[TrilinosImpl::solve] Call trilinos_global_solve_ " << std::endl;
  #endif

  if (consts::trilinos_preconditioners.count(preconditioner_) == 0) {
    auto prec_name = consts::preconditioner_type_to_name.at(preconditioner_); 
    throw std::runtime_error("[TrilinosLinearAlgebra::solve] ERROR: '" + prec_name + "' is not a valid Trilinos preconditioner.");
  }

  trilinos_global_solve_(trilinos_, Val.data(), R.data(), R_.data(), W_.data(),
                         lEq.FSILS.RI.fNorm, lEq.FSILS.RI.iNorm,
                         lEq.FSILS.RI.itr, lEq.FSILS.RI.callD, lEq.FSILS.RI.dB,
                         lEq.FSILS.RI.success, solver_type, lEq.FSILS.RI.relTol,
                         lEq.FSILS.RI.mItr, lEq.FSILS.RI.sD, prec_type);

  for (int a = 0; a < com_mod.tnNo; a++) {
    for (int i = 0; i < com_mod.R.nrows(); i++) {
      com_mod.R(i,a) = R_(i,com_mod.lhs.map(a));
    }
  } 
}

/// @brief Solve a system of linear equations assembled by Trilinos.
void TrilinosLinearAlgebra::TrilinosImpl::solve_assembled(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res)
{
  lEq.FSILS.RI.success = false;
  int solver_type = static_cast<int>(lEq.ls.LS_type);
  int prec_type = static_cast<int>(preconditioner_);
  bool assembled = true;
  #define n_debug_solve_assembled
  #ifdef debug_solve_assembled
  auto prec_name = consts::preconditioner_type_to_name.at(preconditioner_); 
  std::cout << "[TrilinosImpl::solve_assembled] ---------- solve_assembled ---------- " << std::endl;
  std::cout << "[TrilinosImpl::solve_assembled] solver_type: " << solver_type << std::endl;
  std::cout << "[TrilinosImpl::solve_assembled] prec_type: " << prec_name << std::endl;
  std::cout << "[TrilinosImpl::solve_assembled] R_.size(): " << R_.size() << std::endl;
  std::cout << "[TrilinosImpl::solve_assembled] W_.size(): " << W_.size() << std::endl;
  std::cout << "[TrilinosImpl::solve_assembled] lEq.FSILS.RI.success: "
            << lEq.FSILS.RI.success << std::endl;
  #endif

  if (consts::trilinos_preconditioners.count(preconditioner_) == 0) {
    auto prec_name = consts::preconditioner_type_to_name.at(preconditioner_);
    throw std::runtime_error("[TrilinosLinearAlgebra::solve_assembled] ERROR: '" + prec_name + "' is not a valid Trilinos preconditioner.");
  }

  init_dir_and_coup_neu(com_mod, incL, res);

  trilinos_solve_(trilinos_, R_.data(), W_.data(), lEq.FSILS.RI.fNorm,
                  lEq.FSILS.RI.iNorm, lEq.FSILS.RI.itr, lEq.FSILS.RI.callD,
                  lEq.FSILS.RI.dB, lEq.FSILS.RI.success, solver_type,
                  lEq.FSILS.RI.relTol, lEq.FSILS.RI.mItr, lEq.FSILS.RI.sD,
                  prec_type, assembled);

  for (int a = 0; a < com_mod.tnNo; a++) {
    for (int i = 0; i < com_mod.R.nrows(); i++) {
      com_mod.R(i,a) = R_(i,com_mod.lhs.map(a));
    }
  }
}