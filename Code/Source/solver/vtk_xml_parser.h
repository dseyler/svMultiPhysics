// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ComMod.h"

#ifndef VTK_XML_PARSER
#define VTK_XML_PARSER 

namespace vtk_xml_parser {

enum class VtkFileFormat {
  VTP,
  VTU
};

class VtkFileExtentions {
  public:
    const static std::string VTK_VTU_EXTENSION;
    const static std::string VTK_VTP_EXTENSION;
};

void load_fiber_direction_vtu(const std::string& file_name, const std::string& data_name, const int idx,
    const int nsd, mshType& mesh);

/// @brief Which optional groups of per-element active-stress params a spatial VTU supplied.
struct ActiveStressDistInfo {
  bool has_eta = false;     // all three of eta_f/eta_s/eta_n present
  bool has_delay = false;   // delay present
  bool any() const { return has_eta || has_delay; }
};

ActiveStressDistInfo load_active_stress_directional_distribution_vtu(const std::string& file_name,
    Array<double>& elemental_distribution);

void load_vtp(const std::string& file_name, faceType& face);

void load_vtp(const std::string& file_name, mshType& mesh);

void load_vtu(const std::string& file_name, mshType& mesh);

void load_vtu(const std::string& file_name, faceType& face);

void load_time_varying_field_vtu(const std::string file_name, const std::string field_name, mshType& mesh);

};

#endif


