1. Get dependencies.

   ```bash
   cd <mypath>
   wget http://cern.ch/biodynamo-lfs/third-party/metis-5.1.0.tar # METIS
   git clone https://github.com/hypre-space/hypre.git # HYPRE
   git clone https://github.com/mfem/mfem.git #MFEM
   ```

2. Build and install METIS

   ```bash
   cd <mypath>
   export myCC=<mycompiler, e.g., gcc, clang>
   export myCXX=<mycompiler, e.g., g++, clang++>
   tar -xf metis-5.1.0.tar && \
   cd metis-5.1.0 && \
   make config shared=0 cc=$myCC cxx=$myCXX && \
   make && \
   sudo make install  # last one is optional
   ```

3. Build and install HYPRE

   ```bash
   cd <mypath>
   cd hypre/src && \
   mkdir build && \
   cd build && \
   cmake .. && \
   make -j $(nproc --all) && \
   sudo make install # last one is optional
   ```

4. Build and install MFEM

   ```bash
   cd <mypath>
   cd mfem && \
   mkdir build && \
   cd build && \
   cmake -DMFEM_USE_MPI=ON -DMFEM_USE_LAPACK=ON .. && \
   make -j $(nproc --all) && \
   sudo make install # last one is optional
   ```

5. Test MFEM installation

   ```bash 
   cd <mypath>
   cd mfem/build/examples && \
   make ex1p && \
   mpirun -np 4 ./ex1p
   ```
