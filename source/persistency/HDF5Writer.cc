// ----------------------------------------------------------------------------
// nexus | HDF5Writer.cc
//
// This class writes the h5 nexus output file.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "HDF5Writer.h"

#include <sstream>
#include <cstring>
#include <stdlib.h>
#include <vector>

#include <stdint.h>
#include <iostream>

using namespace nexus;


HDF5Writer::HDF5Writer():
  file_(0), irun_(0), ismp_(0), ihit_(0),
  ipart_(0), ipos_(0), istep_(0)
{
}

HDF5Writer::~HDF5Writer()
{
}

void HDF5Writer::Open(std::string fileName, bool debug)
{
  firstEvent_= true;

  file_ = H5Fcreate( fileName.c_str(), H5F_ACC_TRUNC,
                      H5P_DEFAULT, H5P_DEFAULT );

  std::string group_name = "/MC";
  size_t group = createGroup(file_, group_name);

  std::string run_table_name = "configuration";
  memtypeRun_ = createRunType();
  runTable_ = createTable(group, run_table_name, memtypeRun_);

  std::string sns_data_table_name = "sns_response";
  memtypeSnsData_ = createSensorDataType();
  snsDataTable_ = createTable(group, sns_data_table_name, memtypeSnsData_);

  std::string hit_info_table_name = "hits";
  memtypeHitInfo_ = createHitInfoType();
  hitInfoTable_ = createTable(group, hit_info_table_name, memtypeHitInfo_);

  std::string particle_info_table_name = "particles";
  memtypeParticleInfo_ = createParticleInfoType();
  particleInfoTable_ = createTable(group, particle_info_table_name, memtypeParticleInfo_);

  std::string sns_pos_table_name = "sns_positions";
  memtypeSnsPos_ = createSensorPosType();
  snsPosTable_ = createTable(group, sns_pos_table_name, memtypeSnsPos_);

  if (debug) {
    std::string debug_group_name = "/DEBUG";
    size_t debug_group = createGroup(file_, debug_group_name);
    std::string step_table_name = "steps";
    memtypeStep_ = createStepType();
    stepTable_   = createTable(debug_group, step_table_name, memtypeStep_);
  }

  isOpen_ = true;
}

void HDF5Writer::Close()
{
  isOpen_=false;
  H5Fclose(file_);
}

void HDF5Writer::WriteRunInfo(const char* param_key, const char* param_value)
{
  run_info_t runData;
  memset(runData.param_key,   0, CONFLEN);
  memset(runData.param_value, 0, CONFLEN);
  strcpy(runData.param_key, param_key);
  strcpy(runData.param_value, param_value);
  writeRun(&runData, runTable_, memtypeRun_, irun_);

  irun_++;
}


void HDF5Writer::WriteSensorDataInfo(int64_t evt_number, unsigned int sensor_id, unsigned int time_bin, unsigned int charge)
{
  sns_data_t snsData;
  snsData.event_id = evt_number;
  snsData.sensor_id = sensor_id;
  snsData.time_bin = time_bin;
  snsData.charge = charge;
  writeSnsData(&snsData, snsDataTable_, memtypeSnsData_, ismp_);

  ismp_++;
}

void HDF5Writer::WriteHitInfo(int64_t evt_number, int particle_indx, int hit_indx, float hit_position_x, float hit_position_y, float hit_position_z, float hit_time, float hit_energy, const char* label)
{
  hit_info_t trueInfo;
  trueInfo.event_id = evt_number;
  trueInfo.x = hit_position_x;
  trueInfo.y = hit_position_y;
  trueInfo.z = hit_position_z;
  trueInfo.time = hit_time;
  trueInfo.energy = hit_energy;
  memset(trueInfo.label, 0, STRLEN);
  strcpy(trueInfo.label, label);
  trueInfo.particle_id = particle_indx;
  trueInfo.hit_id = hit_indx;
  writeHit(&trueInfo,  hitInfoTable_, memtypeHitInfo_, ihit_);

  ihit_++;
}

void HDF5Writer::WriteParticleInfo(int64_t evt_number, int particle_indx, const char* particle_name, char primary, int mother_id, float initial_vertex_x, float initial_vertex_y, float initial_vertex_z, float initial_vertex_t, float final_vertex_x, float final_vertex_y, float final_vertex_z, float final_vertex_t, const char* initial_volume, const char* final_volume, float ini_momentum_x, float ini_momentum_y, float ini_momentum_z, float final_momentum_x, float final_momentum_y, float final_momentum_z, float kin_energy, float length, const char* creator_proc, const char* final_proc)
{
  particle_info_t trueInfo;
  trueInfo.event_id = evt_number;
  trueInfo.particle_id = particle_indx;
  memset(trueInfo.particle_name, 0, STRLEN);
  strcpy(trueInfo.particle_name, particle_name);
  trueInfo.primary = primary;
  trueInfo.mother_id = mother_id;
  trueInfo.initial_x = initial_vertex_x;
  trueInfo.initial_y = initial_vertex_y;
  trueInfo.initial_z = initial_vertex_z;
  trueInfo.initial_t = initial_vertex_t;
  trueInfo.final_x = final_vertex_x;
  trueInfo.final_y = final_vertex_y;
  trueInfo.final_z = final_vertex_z;
  trueInfo.final_t = final_vertex_t;
  memset(trueInfo.initial_volume, 0, STRLEN);
  strcpy(trueInfo.initial_volume, initial_volume);
  memset(trueInfo.final_volume, 0, STRLEN);
  strcpy(trueInfo.final_volume, final_volume);
  trueInfo.initial_momentum_x = ini_momentum_x;
  trueInfo.initial_momentum_y = ini_momentum_y;
  trueInfo.initial_momentum_z = ini_momentum_z;
  trueInfo.final_momentum_x = final_momentum_x;
  trueInfo.final_momentum_y = final_momentum_y;
  trueInfo.final_momentum_z = final_momentum_z;
  trueInfo.kin_energy = kin_energy;
  trueInfo.length = length;
  memset(trueInfo.creator_proc, 0, STRLEN);
  strcpy(trueInfo.creator_proc, creator_proc);
  memset(trueInfo.final_proc, 0, STRLEN);
  strcpy(trueInfo.final_proc, final_proc);
  writeParticle(&trueInfo,  particleInfoTable_, memtypeParticleInfo_, ipart_);

  ipart_++;
}

void HDF5Writer::WriteSensorPosInfo(unsigned int sensor_id, const char* sensor_name, float x, float y, float z)
{
  sns_pos_t snsPos;
  snsPos.sensor_id = sensor_id;
  memset(snsPos.sensor_name, 0, STRLEN);
  strcpy(snsPos.sensor_name, sensor_name);
  snsPos.x = x;
  snsPos.y = y;
  snsPos.z = z;
  writeSnsPos(&snsPos, snsPosTable_, memtypeSnsPos_, ipos_);

  ipos_++;
}

void HDF5Writer::WriteStep(int64_t evt_number,
                           int particle_id, const char* particle_name,
                           int step_id,
                           const char* initial_volume,
                           const char*   final_volume,
                           const char*      proc_name,
                           float initial_x, float initial_y, float initial_z,
                           float   final_x, float   final_y, float   final_z,
                           float time)
{
  step_info_t step;
  step.event_id    = evt_number;
  step.particle_id = particle_id;
  memset(step.particle_name , 0,  STRLEN);
  strcpy(step.particle_name ,  particle_name);
  step.step_id    = step_id;
  memset(step.initial_volume, 0, STRLEN);
  strcpy(step.initial_volume, initial_volume);
  memset(step.  final_volume, 0, STRLEN);
  strcpy(step.  final_volume,   final_volume);
  memset(step.     proc_name, 0, STRLEN);
  strcpy(step.     proc_name,      proc_name);
  step.initial_x   = initial_x;
  step.initial_y   = initial_y;
  step.initial_z   = initial_z;
  step.  final_x   =   final_x;
  step.  final_y   =   final_y;
  step.  final_z   =   final_z;
  step.time        =      time;

  writeStep(&step, stepTable_, memtypeStep_, istep_);

  istep_++;
}
