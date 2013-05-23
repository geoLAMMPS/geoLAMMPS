# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp -p compute_coord_gran.cpp ..
  cp -p fix_fluiddrag.cpp ..
  cp -p fix_freeze.cpp ..
  cp -p fix_pour.cpp ..
  cp -p fix_wall_gran.cpp ..
  cp -p fix_write_insurance_shear_history.cpp ..
  cp -p pair_gran_hertz_history.cpp ..
  cp -p pair_gran_hooke.cpp ..
  cp -p pair_gran_hooke_history.cpp ..
  cp -p pair_gran_shm_history.cpp ..

  cp -p compute_coord_gran.h ..
  cp -p fix_fluiddrag.h ..
  cp -p fix_freeze.h ..
  cp -p fix_pour.h ..
  cp -p fix_wall_gran.h ..
  cp -p fix_write_insurance_shear_history.h ..
  cp -p pair_gran_hertz_history.h ..
  cp -p pair_gran_hooke.h ..
  cp -p pair_gran_hooke_history.h ..
  cp -p pair_gran_shm_history.h ..

elif (test $1 = 0) then

  rm -f ../compute_coord_gran.cpp
  rm -f ../fix_fluiddrag.cpp
  rm -f ../fix_freeze.cpp
  rm -f ../fix_pour.cpp
  rm -f ../fix_wall_gran.cpp
  rm -f ../fix_write_insurance_shear_history.cpp
  rm -f ../pair_gran_hertz_history.cpp
  rm -f ../pair_gran_hooke.cpp
  rm -f ../pair_gran_hooke_history.cpp
  rm -f ../pair_gran_shm_history.cpp

  rm -f ../compute_coord_gran.h
  rm -f ../fix_fluiddrag.h
  rm -f ../fix_freeze.h
  rm -f ../fix_pour.h
  rm -f ../fix_wall_gran.h
  rm -f ../fix_write_insurance_shear_history.h
  rm -f ../pair_gran_hertz_history.h
  rm -f ../pair_gran_hooke.h
  rm -f ../pair_gran_hooke_history.h
  rm -f ../pair_gran_shm_history.h

fi
