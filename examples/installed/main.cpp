// SPDX-License-Identifier: GPL-3.0-or-later
#include <cleanseams/build_info.hpp>
#include <cleanseams/mesh_analysis.hpp>
#include <iostream>
int main() {
    cleanseams::MeshInput mesh;
    mesh.positions={{0,0,0},{2,0,0},{2,2,0},{0,2,0}};
    mesh.face_offsets={0,4};mesh.corner_vertices={0,1,2,3};
    const auto analysis=cleanseams::analyze_mesh(mesh);
    if(analysis.triangle_corners.size()!=2 || analysis.face_area_3d.size()!=1 ||
       analysis.face_area_3d[0]!=4 || analysis.corner_faces.size()!=4) return 1;
    const auto info=cleanseams::build_info();
    if(!info.at("features").at("mesh_analysis") || info.at("features").at("python_engine")) return 2;
    std::cout<<info.dump()<<'\n';
}
