//----------------------------------------------------------------------------------------------
// This file is part of ESEPP (version 1.5).
//
// ESEPP is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ESEPP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ESEPP.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright (c) Alexander Gramolin, 2019. E-mail: gramolin (at) bu.edu
// https://github.com/gramolin/esepp/
//==============================================================================================


int Dialog ()
{
char mychar[64];

//----------------------------------------------------------------------------------------------
// Welcome message:
std::cout << std::endl;
std::cout << "/////////////////////////////////////////////////////////////////////////////////" << std::endl;
std::cout << "//    ESEPP v1.5:  Elastic Scattering of Electrons and Positrons on Protons    //" << std::endl;
std::cout << "//     Copyright (c) Alexander Gramolin, 2019.     E-mail: gramolin"; std::cout << "@bu.edu     //" << std::endl;
std::cout << "//     https://github.com/gramolin/esepp/   http://arxiv.org/abs/1401.2959     //" << std::endl;
std::cout << "//                                                                             //" << std::endl;
std::cout << "// This program comes with absolutely no warranty; for details run \"esepp -w\". //" << std::endl;
std::cout << "//                                                                             //" << std::endl;
std::cout << "// This is free software, and you are welcome to redistribute it under certain //" << std::endl;
std::cout << "//     conditions; run the program with a key -c (\"esepp -c\") for details.     //" << std::endl;
std::cout << "/////////////////////////////////////////////////////////////////////////////////" << std::endl;
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Type of scattered leptons:
if (flag_vepp == false)
  {
  std::cout << std::endl;
  std::cout << "*PLEASE SELECT THE TYPE OF SCATTERED LEPTONS:" << std::endl;
  std::cout << " 1 -- electrons (e-) only;" << std::endl;
  std::cout << " 2 -- positrons (e+) only;" << std::endl;
  std::cout << " 3 -- both electrons and positrons (e- and e+, default);" << std::endl << std::endl;
  std::cout << " 4 -- muons (mu-) only;" << std::endl;
  std::cout << " 5 -- antimuons (mu+) only;" << std::endl;
  std::cout << " 6 -- both muons and antimuons (mu- and mu+)." << std::endl;
  std::cout << std::endl << " Your choice (3): " << std::flush;

  while (1)
    {
    mychar[0] = 0; mychar[1] = 0;
    std::cin.getline(mychar, 64);
    if (mychar[1] != 0) { std::cout << " Wrong input! Try again: " << std::flush; continue; }
    switch (mychar[0])
      {
      case '1': { flag_lepton = 1; m = m_e; break; }
      case '2': { flag_lepton = 2; m = m_e; break; }
      case '3': case 0: { flag_lepton = 3; m = m_e; break; }
      case '4': { flag_lepton = 4; m = m_mu; break; }
      case '5': { flag_lepton = 5; m = m_mu; break; }
      case '6': { flag_lepton = 6; m = m_mu; break; }
      default: { std::cout << " Wrong input! Try again: " << std::flush; continue; }
      }
    m2 = Pow2(m); m4 = Pow4(m);
    break;
    }
  }
else { flag_lepton = 3; m = m_e; }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Events in accordance with the Rosenbluth formula
std::cout << std::endl;
std::cout << "*DO YOU NEED EVENTS IN ACCORDANCE WITH THE ROSENBLUTH FORMULA?" << std::endl;
std::cout << " 1 -- yes;" << std::endl;
std::cout << " 2 -- no (default)." << std::endl;
std::cout << std::endl << " Your choice (2): " << std::flush;

while (1)
  {
  mychar[0] = 0; mychar[1] = 0;
  std::cin.getline(mychar, 64);
  if (mychar[1] != 0) { std::cout << " Wrong input! Try again: " << std::flush; continue; }
  switch (mychar[0])
    {
    case '1': { flag_rosen = true; break; }
    case '2': case 0: { flag_rosen = false; break; }
    default: { std::cout << " Wrong input! Try again: " << std::flush; continue; }
    }
  break;
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Internal structure of the proton:
std::cout << std::endl;
std::cout << "*INTERNAL STRUCTURE OF THE PROTON:" << std::endl;
std::cout << " 1 -- a point-like proton;" << std::endl;
std::cout << " 2 -- a proton with the dipole form factors (default);" << std::endl;
std::cout << " 3 -- a proton with the Kelly form factors;" << std::endl;
std::cout << " 4 -- a proton with the Puckett form factors;" << std::endl;
std::cout << " 5 -- specify your own parametrization (see the file \"const.h\")." << std::endl;
std::cout << std::endl << " Your choice (2): " << std::flush;

while (1)
  {
  mychar[0] = 0; mychar[1] = 0;
  std::cin.getline(mychar, 64);
  if (mychar[1] != 0) { std::cout << " Wrong input! Try again: " << std::flush; continue; }
  switch (mychar[0])
    {
    case '1': { flag_struct = 1; break; }
    case '2': case 0: { flag_struct = 2; break; }
    case '3': { flag_struct = 3; break; }
    case '4': { flag_struct = 4; break; }
    case '5': { flag_struct = 5; break; }
    default: { std::cout << " Wrong input! Try again: " << std::flush; continue; }
    }
  break;
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Mode of calculation for bremsstrahlung:
std::cout << std::endl;
std::cout << "*PLEASE SELECT ONE OF THE MODES FOR BREMSSTRAHLUNG:" << std::endl;
std::cout << " 1 -- primary soft-photon approximation, only lepton bremsstrahlung;" << std::endl;
std::cout << " 2 -- primary soft-photon approximation, only proton bremsstrahlung;" << std::endl;
std::cout << " 3 -- primary soft-photon approximation, all the terms;" << std::endl << std::endl;
std::cout << " 4 -- modified soft-photon approximation, only lepton bremsstrahlung;" << std::endl;
std::cout << " 5 -- modified soft-photon approximation, only proton bremsstrahlung;" << std::endl;
std::cout << " 6 -- modified soft-photon approximation, all the terms;" << std::endl << std::endl;
std::cout << " 7 -- improved soft-photon approximation, only lepton bremsstrahlung;" << std::endl;
std::cout << " 8 -- improved soft-photon approximation, only proton bremsstrahlung;" << std::endl;
std::cout << " 9 -- improved soft-photon approximation, all the terms;" << std::endl << std::endl;
std::cout << "10 -- accurate QED calculation, only lepton bremsstrahlung;" << std::endl;
std::cout << "11 -- accurate QED calculation, only proton bremsstrahlung;" << std::endl;
std::cout << "12 -- accurate QED calculation, all the terms (default)." << std::endl;
std::cout << std::endl << " Your choice (12): " << std::flush;

while (1)
  {
  mychar[0] = 0; mychar[1] = 0;
  std::cin.getline(mychar, 64);
  if (mychar[0] == 0) { flag_mode = 12; break; } // Default
  if (atoi(mychar) > 0 && atoi(mychar) < 13) { flag_mode = atoi(mychar); break; }
  else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Vacuum polarization correction:
std::cout << std::endl;
std::cout << "*PLEASE SELECT ONE OF THE MODES FOR THE VACUUM POLARIZATION:" << std::endl;
std::cout << " 1 -- only electron-positron loops;" << std::endl;
std::cout << " 2 -- full leptonic contribution;" << std::endl;
std::cout << " 3 -- full vacuum polarization correction (default)." << std::endl;
std::cout << std::endl << " Your choice (3): " << std::flush;
  
while (1)
  {
  mychar[0] = 0; mychar[1] = 0;
  std::cin.getline(mychar, 64);
  if (mychar[1] != 0) { std::cout << " Wrong input! Try again: " << std::flush; continue; }
  switch (mychar[0])
    {
    case '1': { flag_vpol = 1; break; }
    case '2': { flag_vpol = 2; break; }
    case '3': case 0:
      {
      flag_vpol = 3;
      
      fvpol = fopen("vpol.dat", "r"); // Opening the file "vpol.dat"
      if (fvpol == NULL) { std::cout << " Can't open file \"vpol.dat\"!" << std::endl; return EXIT_FAILURE; }
      
      break;
      }
    default: { std::cout << " Wrong input! Try again: " << std::flush; continue; }
    }
  break;
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Two-photon exchange amplitudes:
std::cout << std::endl;
std::cout << "*PLEASE SELECT ONE OF THE MODES FOR THE TWO-PHOTON EXCHANGE AMPLITUDES:" << std::endl;
std::cout << " 1 -- approach of Mo & Tsai (default);" << std::endl;
std::cout << " 2 -- approach of Maximon & Tjon." << std::endl;
std::cout << std::endl << " Your choice (1): " << std::flush;
  
while (1)
  {
  mychar[0] = 0; mychar[1] = 0;
  std::cin.getline(mychar, 64);
  if (mychar[1] != 0) { std::cout << " Wrong input! Try again: " << std::flush; continue; }
  switch (mychar[0])
    {
    case '1': case 0: { flag_tpe = 1; break; }
    case '2': { flag_tpe = 2; break; }
    default: { std::cout << " Wrong input! Try again: " << std::flush; continue; }
    }
  break;
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Kinematic parameters:
std::cout << std::endl;
std::cout << "*KINEMATIC PARAMETERS:" << std::endl;
std::cout << " Full energy of incident leptons, MeV (1000): " << std::flush;

while (1)
  {
  mychar[0] = 0; mychar[1] = 0;
  std::cin.getline(mychar, 64);
  if (mychar[0] == 0) { E_li = 1.; break; } // Default: E_li = 1 GeV
  else
    {
    E_li = atof(mychar)/1000.;
    if (E_li > m) break;
    else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
    }
  }

if (E_li > 0.05) std::cout << " Cut-off energy for bremsstrahlung photons, MeV (1): " << std::flush;
else std::cout << " Cut-off energy for bremsstrahlung photons, MeV (0.1): " << std::flush;

while (1)
  {
  mychar[0] = 0; mychar[1] = 0;
  std::cin.getline(mychar, 64);
  if (mychar[0] == 0) { E_g_cut = 0.0001; if (E_li > 0.05) E_g_cut = 0.001; break; } // Default
  else
    {
    E_g_cut = atof(mychar)/1000.;
    if (1000.*E_g_cut >= 0.01 && E_g_cut < E_li - m) break; // E_g_cut must be not less than 0.01 MeV
    else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
    }
  }

std::cout << " Maximum energy for bremsstrahlung photons, MeV (" << int (700*(E_li - m) + 0.5) << "): " << std::flush;

while (1)
  {
  mychar[0] = 0; mychar[1] = 0;
  std::cin.getline(mychar, 64);
  if (mychar[0] == 0) { E_g_max = (int (700*(E_li - m) + 0.5))/1000.; break; } // Default: E_g_max = 70% of (E_li - m)
  else
    {
    E_g_max = atof(mychar)/1000.;
    if (E_g_max > E_g_cut && E_g_max <= E_li - m) break;
    else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
    }
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Polar angles of the scattered lepton (theta):
std::cout << std::endl;
std::cout << "*POLAR ANGLES OF THE SCATTERED LEPTON (THETA) IN THE LAB FRAME:" << std::endl;
std::cout << " 1 -- Full angular range (from 1 to 180 degrees);" << std::endl;
if (flag_vepp == false) std::cout << " 2 -- To specify any angular range (default)." << std::endl;
if (flag_vepp == true)
  {
  std::cout << " 2 -- To specify any angular range (default);" << std::endl << std::endl;
  std::cout << " 3 -- VEPP-3, Run I, SA (from 5 to 25 degrees);" << std::endl;
  std::cout << " 4 -- VEPP-3, Run I, MA (from 13 to 40 degrees);" << std::endl;
  std::cout << " 5 -- VEPP-3, Run I, LA (from 45 to 95 degrees);" << std::endl << std::endl;
  std::cout << " 6 -- VEPP-3, Run II, MA (from 13 to 40 degrees);" << std::endl;
  std::cout << " 7 -- VEPP-3, Run II, LA (from 60 to 115 degrees);" << std::endl << std::endl;
  std::cout << " 8 -- VEPP-3, Run III, MA (from 17 to 40 degrees);" << std::endl;
  std::cout << " 9 -- VEPP-3, Run III, LA (from 60 to 115 degrees)." << std::endl;
  }
std::cout << std::endl << " Your choice (2): " << std::flush;

while (1)
  {
  mychar[0] = 0; mychar[1] = 0;
  std::cin.getline(mychar, 64);
  if (mychar[1] != 0) { std::cout << " Wrong input! Try again: " << std::flush; continue; }
  switch (mychar[0])
    {
    case '1': { theta_min = 1.*degrad; theta_max = Pi; break; }
    case '2': case 0:
      {
      std::cout << std::endl << " Minimum value, degree (1): " << std::flush;
      while (1)
        {
        mychar[0] = 0; mychar[1] = 0;
        std::cin.getline(mychar, 64);
        if (mychar[0] == 0) { theta_min = 1.*degrad; break; }
        else
          {
          theta_min = atof(mychar);
          if (theta_min > 0 && theta_min < 180) { theta_min *= degrad; break; }
          else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
          }
        }

      std::cout << " Maximum value, degree (180): " << std::flush;
      while (1)
        {
        mychar[0] = 0; mychar[1] = 0;
        std::cin.getline(mychar, 64);
        if (mychar[0] == 0) { theta_max = Pi; break; }
        else
          {
          theta_max = atof(mychar);
          if (theta_max*degrad > theta_min && theta_max <= 180) { theta_max *= degrad; break; }
          else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
          }
        }

      break;
      }
    case '3': { if (flag_vepp == false) { std::cout << " Wrong input! Try again: " << std::flush; continue; } theta_min = 5.*degrad; theta_max = 25.*degrad; break; }
    case '4': { if (flag_vepp == false) { std::cout << " Wrong input! Try again: " << std::flush; continue; } theta_min = 13.*degrad; theta_max = 40.*degrad; break; }
    case '5': { if (flag_vepp == false) { std::cout << " Wrong input! Try again: " << std::flush; continue; } theta_min = 45.*degrad; theta_max = 95.*degrad; break; }
    case '6': { if (flag_vepp == false) { std::cout << " Wrong input! Try again: " << std::flush; continue; } theta_min = 13.*degrad; theta_max = 40.*degrad; break; }
    case '7': { if (flag_vepp == false) { std::cout << " Wrong input! Try again: " << std::flush; continue; } theta_min = 60.*degrad; theta_max = 115.*degrad; break; }
    case '8': { if (flag_vepp == false) { std::cout << " Wrong input! Try again: " << std::flush; continue; } theta_min = 17.*degrad; theta_max = 40.*degrad; break; }
    case '9': { if (flag_vepp == false) { std::cout << " Wrong input! Try again: " << std::flush; continue; } theta_min = 60.*degrad; theta_max = 115.*degrad; break; }
    default: { std::cout << " Wrong input! Try again: " << std::flush; continue; }
    }
  break;
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Azimuthal angles of the scattered lepton (phi):
std::cout << std::endl;
std::cout << "*AZIMUTHAL ANGLES OF THE SCATTERED LEPTON (PHI):" << std::endl;
std::cout << " 1 -- convention where 0 < PHI < 360 degrees;" << std::endl;
std::cout << " 2 -- convention where -180 < PHI < 180 degrees (default)." << std::endl;
std::cout << std::endl << " Your choice (2): " << std::flush;

while (1)
  {
  mychar[0] = 0; mychar[1] = 0;
  std::cin.getline(mychar, 64);
  if (mychar[1] != 0) { std::cout << " Wrong input! Try again: " << std::flush; continue; }
  switch (mychar[0])
    {
    case '1':
      {
      flag_phi = 1;

      std::cout << std::endl << " Minimum value, degree (0): " << std::flush;
      while (1)
        {
        mychar[0] = 0; mychar[1] = 0;
        std::cin.getline(mychar, 64);
        if (mychar[0] == 0) { phi_min = 0.; break; }
        else
          {
          phi_min = atof(mychar);
          if ((phi_min >= 0 || mychar[0] == '0') && phi_min < 360) { phi_min *= degrad; break; }
          else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
          }
        }

      if (flag_vepp == false) std::cout << " Maximum value, degree (360): " << std::flush;
      else std::cout << " Maximum value, degree (70): " << std::flush;
      while (1)
        {
        mychar[0] = 0; mychar[1] = 0;
        std::cin.getline(mychar, 64);
        if (mychar[0] == 0)
	  {
	  if (flag_vepp == false) { phi_max = 360.*degrad; break; }
	  else { phi_max = 70.*degrad; break; }
	  }
        else
          {
          phi_max = atof(mychar);
          if (phi_max*degrad > phi_min && phi_max <= 360) { phi_max *= degrad; break; }
          else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
          }
        }

      break;
      }

    case '2': case 0:
      {
      flag_phi = 2;

      if (flag_vepp == false) std::cout << std::endl << " Minimum value, degree (-180): " << std::flush;
      else std::cout << std::endl << " Minimum value, degree (-35): " << std::flush;
      while (1)
        {
        mychar[0] = 0; mychar[1] = 0;
        std::cin.getline(mychar, 64);
        if (mychar[0] == 0)
	  {
	  if (flag_vepp == false) { phi_min = -180.*degrad; break; }
	  else  { phi_min = -35.*degrad; break; }
	  }
        else
          {
          phi_min = atof(mychar);
          if ((phi_min != 0 || mychar[0] == '0') && phi_min >= -180 && phi_min < 180) { phi_min *= degrad; break; }
          else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
          }
        }

      if (flag_vepp == false) std::cout << " Maximum value, degree (180): " << std::flush;
      else std::cout << " Maximum value, degree (35): " << std::flush;
      while (1)
        {
        mychar[0] = 0; mychar[1] = 0;
        std::cin.getline(mychar, 64);
        if (mychar[0] == 0)
	  {
	  if (flag_vepp == false) { phi_max = 180.*degrad; break; }
	  else { phi_max = 35.*degrad; break; }
	  }
        else
          {
          phi_max = atof(mychar);
          if (phi_max*degrad > phi_min && phi_max <= 180) { phi_max *= degrad; break; }
          else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
          }
        }

      break;
      }
    default: { std::cout << " Wrong input! Try again: " << std::flush; continue; }
    }
  break;
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Storage cell parameters:
if (flag_target == true)
  {
  std::cout << std::endl;
  std::cout << "*STORAGE CELL PARAMETERS:" << std::endl;
  std::cout << " 1 -- uniform distribution of gas pressure;" << std::endl;
  std::cout << " 2 -- triangle distribution of gas pressure (default)." << std::endl;
  std::cout << std::endl << " Your choice (2): " << std::flush;

  while (1)
    {
    mychar[0] = 0; mychar[1] = 0;
    std::cin.getline(mychar, 64);
    if (mychar[1] != 0) { std::cout << " Wrong input! Try again: " << std::flush; continue; }
    switch (mychar[0])
      {
      case '1': { flag_cell = 1; break; }
      case '2': case 0: { flag_cell = 2; break; }
      default: { std::cout << " Wrong input! Try again: " << std::flush; continue; }
      }
    break;
    }

  std::cout << std::endl << " Minimum Z coordinate of the storage cell, mm (-200): " << std::flush;

  while (1)
    {
    mychar[0] = 0; mychar[1] = 0;
    std::cin.getline(mychar, 64);
    if (mychar[0] == 0) { cell_zmin = -200.; break; } // Default: cell_zmin = -200 mm
    else
      {
      cell_zmin = atof(mychar);
      if (cell_zmin != 0 || mychar[0] == '0') break;
      else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
      }
    }

  std::cout << " Maximum Z coordinate of the storage cell, mm (200): " << std::flush;

  while (1)
    {
    mychar[0] = 0; mychar[1] = 0;
    std::cin.getline(mychar, 64);
    if (mychar[0] == 0) { cell_zmax = 200.; break; } // Default: cell_zmax = 200 mm
    else
      {
      cell_zmax = atof(mychar);
      if ((cell_zmax != 0 || mychar[0] == '0') && cell_zmax > cell_zmin) break;
      else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
      }
    }

  if (flag_cell == 2)
    {
    std::cout << " Z coordinate of the point of gas injection, mm (0): " << std::flush;

    while (1)
      {
      mychar[0] = 0; mychar[1] = 0;
      std::cin.getline(mychar, 64);
      if (mychar[0] == 0) { cell_zinj = 0.; break; } // Default: cell_zinj = 0 mm
      else
        {
        cell_zinj = atof(mychar);
        if ((cell_zinj != 0 || mychar[0] == '0') && cell_zinj >= cell_zmin && cell_zinj <= cell_zmax) break;
        else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
        }
      }
    }
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// To specify the beam current integral instead of the number of events (for VEPP-3 only):
if (flag_vepp == true)
  {
  std::cout << std::endl;
  std::cout << "*DO YOU WANT TO SPECIFY THE BEAM CURRENT INTEGRAL REQUIRED?" << std::endl;
  std::cout << " 1 -- yes (default);" << std::endl;
  std::cout << " 2 -- no." << std::endl;
  std::cout << std::endl << " Your choice (1): " << std::flush;

  while (1)
    {
    mychar[0] = 0; mychar[1] = 0;
    std::cin.getline(mychar, 64);
    if (mychar[1] != 0) { std::cout << " Wrong input! Try again: " << std::flush; continue; }
    switch (mychar[0])
      {
      case '1': case 0:
        {
	flag_bint = true;
	std::cout << std::endl << " Beam current integral required, kC (10): " << std::flush;
	
        while (1)
          {
          mychar[0] = 0; mychar[1] = 0;
          std::cin.getline(mychar, 64);
          if (mychar[0] == 0) { kiloc = 10.; break; } // Default: kiloc = 10 kC
          else
            {
            kiloc = atof(mychar);
            if (kiloc > 0. && kiloc <= 1000.) break;
            else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
            }
          }
          
        break;
	}
      case '2': { flag_bint = false; break; }
      default: { std::cout << " Wrong input! Try again: " << std::flush; continue; }
      }
    break;
    }
  }
else flag_bint = false;
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Number of events requried:
if (flag_bint == false)
  {
  std::cout << std::endl;
  std::cout << "*PLEASE ENTER THE NUMBER OF EVENTS REQUIRED (1000000): " << std::flush;

  while (1)
    {
    mychar[0] = 0; mychar[1] = 0;
    std::cin.getline(mychar, 64);
    if (mychar[0] == 0) { nevents = 1000000; break; }
    else
      {
      nevents = atol(mychar);
      if (nevents >= 500 && nevents <= 1e9) break;
      else { std::cout << " Wrong input! Try again: " << std::flush; continue; }
      }
    }
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Format of the output files:
std::cout << std::endl;
std::cout << "*FORMAT OF THE OUTPUT FILES:" << std::endl;
std::cout << " 1 -- *.dat files;" << std::endl;
std::cout << " 2 -- *.root files (default);" << std::endl;
std::cout << " 3 -- both *.dat and *.root files." << std::endl;
std::cout << std::endl << " Your choice (2): " << std::flush;

while (1)
  {
  mychar[0] = 0; mychar[1] = 0;
  std::cin.getline(mychar, 64);
  if (mychar[1] != 0) { std::cout << " Wrong input! Try again: " << std::flush; continue; }
  switch (mychar[0])
    {
    case '1': { flag_dat = true; flag_root = false; break; }
    case '2': case 0: { flag_dat = false; flag_root = true; break; }
    case '3': { flag_dat = true; flag_root = true; break; }
    default: { std::cout << " Wrong input! Try again: " << std::flush; continue; }
    }
  break;
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Name (prefix) for the output files:
std::cout << std::endl;
std::cout << "*Please enter the name (prefix) for the output files: " << std::flush;

while (1)
  {
  mychar[0] = 0; mychar[1] = 0;
  std::cin.getline(mychar, 64);
  filename = mychar;

  if (filename == "") filename = "events"; // Default

  if (flag_dat == true) // To produce *.dat files
    {
    // Rosenbluth events:
    if (flag_rosen == true)
      {
      if (flag_lepton < 4) f0 = fopen(filename + "_e0.dat", "w"); // For e-/e+
      else f0 = fopen(filename + "_mu0.dat", "w"); // For mu-/mu+
  
      if (f0 == NULL)
        {
        if (flag_lepton < 4) std::cout << "Error opening file " << filename << "_e0.dat!" << std::endl;
        else std::cout << "Error opening file " << filename << "_mu0.dat!" << std::endl;
        std::cout << "Try again: " << std::flush;
        continue;
        }
      }
      
    // Negatively charged leptons:
    if (flag_lepton == 1 || flag_lepton == 3) fe = fopen(filename + "_e-.dat", "w");
    if (flag_lepton == 4 || flag_lepton == 6) fe = fopen(filename + "_mu-.dat", "w");
    if (flag_lepton != 2 && flag_lepton != 5 && fe == NULL)
      {
      if (flag_lepton == 1 || flag_lepton == 3) std::cout << "Error opening file " << filename << "_e-.dat!" << std::endl;
      if (flag_lepton == 4 || flag_lepton == 6) std::cout << "Error opening file " << filename << "_mu-.dat!" << std::endl;
      std::cout << "Try again: " << std::flush;
      continue;
      }

    // Positively charged leptons:
    if (flag_lepton == 2 || flag_lepton == 3) fp = fopen(filename + "_e+.dat", "w");
    if (flag_lepton == 5 || flag_lepton == 6) fp = fopen(filename + "_mu+.dat", "w");
    if (flag_lepton != 1 && flag_lepton != 4 && fp == NULL)
      {
      if (flag_lepton == 2 || flag_lepton == 3) std::cout << "Error opening file " << filename << "_e+.dat!" << std::endl;
      if (flag_lepton == 5 || flag_lepton == 6) std::cout << "Error opening file " << filename << "_mu+.dat!" << std::endl;
      std::cout << "Try again: " << std::flush;
      continue;
      }
    }

  if (flag_root == true) // To produce *.root files
    {
    // Rosenbluth events:
    if (flag_rosen == true)
      {
      if (flag_lepton < 4) froot_0 = new TFile(filename + "_e0.root", "RECREATE");
      else froot_0 = new TFile(filename + "_mu0.root", "RECREATE");
      if (froot_0->IsZombie())
        {
        if (flag_lepton < 4) std::cout << "Error opening file " << filename << "_e0.root!" << std::endl;
        else std::cout << "Error opening file " << filename << "_mu0.root!" << std::endl;
        std::cout << "Try again: " << std::flush;
        continue;
        }
      if (flag_target == false)
        {
        if (flag_lepton < 4) ntp_0 = new TNtuple("ntp", "Ntuple for e0", "E_l:theta_l:phi_l:E_p:theta_p:phi_p:E_g:theta_g:phi_g");
        else ntp_0 = new TNtuple("ntp", "Ntuple for mu0", "E_l:theta_l:phi_l:E_p:theta_p:phi_p:E_g:theta_g:phi_g");
        }
      else
        {
        if (flag_lepton < 4) ntp_0 = new TNtuple("ntp", "Ntuple for e0", "E_l:theta_l:phi_l:E_p:theta_p:phi_p:E_g:theta_g:phi_g:z");
        else ntp_0 = new TNtuple("ntp", "Ntuple for mu0", "E_l:theta_l:phi_l:E_p:theta_p:phi_p:E_g:theta_g:phi_g:z");
        }
      }
      
    // Negatively charged leptons:
    if (flag_lepton == 1 || flag_lepton == 3) froot_e = new TFile(filename + "_e-.root", "RECREATE");
    if (flag_lepton == 4 || flag_lepton == 6) froot_e = new TFile(filename + "_mu-.root", "RECREATE");
    if (flag_lepton != 2 && flag_lepton != 5 && froot_e->IsZombie())
      {
      if (flag_lepton == 1 || flag_lepton == 3) std::cout << "Error opening file " << filename << "_e-.root!" << std::endl;
      if (flag_lepton == 4 || flag_lepton == 6) std::cout << "Error opening file " << filename << "_mu-.root!" << std::endl;
      std::cout << "Try again: " << std::flush;
      continue;
      }
    if (flag_target == false)
      {
      if (flag_lepton == 1 || flag_lepton == 3) ntp_e = new TNtuple("ntp", "Ntuple for e-", "E_l:theta_l:phi_l:E_p:theta_p:phi_p:E_g:theta_g:phi_g");
      if (flag_lepton == 4 || flag_lepton == 6) ntp_e = new TNtuple("ntp", "Ntuple for mu-", "E_l:theta_l:phi_l:E_p:theta_p:phi_p:E_g:theta_g:phi_g");
      }
    else
      {
      if (flag_lepton == 1 || flag_lepton == 3) ntp_e = new TNtuple("ntp", "Ntuple for e-", "E_l:theta_l:phi_l:E_p:theta_p:phi_p:E_g:theta_g:phi_g:z");
      if (flag_lepton == 4 || flag_lepton == 6) ntp_e = new TNtuple("ntp", "Ntuple for mu-", "E_l:theta_l:phi_l:E_p:theta_p:phi_p:E_g:theta_g:phi_g:z");
      }

    // Positively charged leptons:
    if (flag_lepton == 2 || flag_lepton == 3) froot_p = new TFile(filename + "_e+.root", "RECREATE");
    if (flag_lepton == 5 || flag_lepton == 6) froot_p = new TFile(filename + "_mu+.root", "RECREATE");
    if (flag_lepton != 1 && flag_lepton != 4 && froot_p->IsZombie())
      {
      if (flag_lepton == 2 || flag_lepton == 3) std::cout << "Error opening file " << filename << "_e+.root!" << std::endl;
      if (flag_lepton == 5 || flag_lepton == 6) std::cout << "Error opening file " << filename << "_mu+.root!" << std::endl;
      std::cout << "Try again: " << std::flush;
      continue;
      }
    if (flag_target == false)
      {
      if (flag_lepton == 2 || flag_lepton == 3) ntp_p = new TNtuple("ntp", "Ntuple for e+", "E_l:theta_l:phi_l:E_p:theta_p:phi_p:E_g:theta_g:phi_g");
      if (flag_lepton == 5 || flag_lepton == 6) ntp_p = new TNtuple("ntp", "Ntuple for mu+", "E_l:theta_l:phi_l:E_p:theta_p:phi_p:E_g:theta_g:phi_g");
      }
    else
      {
      if (flag_lepton == 2 || flag_lepton == 3) ntp_p = new TNtuple("ntp", "Ntuple for e+", "E_l:theta_l:phi_l:E_p:theta_p:phi_p:E_g:theta_g:phi_g:z");
      if (flag_lepton == 5 || flag_lepton == 6) ntp_p = new TNtuple("ntp", "Ntuple for mu+", "E_l:theta_l:phi_l:E_p:theta_p:phi_p:E_g:theta_g:phi_g:z");
      }
    }
  break;
  }
//==============================================================================================


std::cout << std::endl << "Data is entered. PLEASE WAIT!" << std::endl;
std::cout << "================================================================================" << std::endl;

return EXIT_SUCCESS;
}

