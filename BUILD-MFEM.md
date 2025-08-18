See also [Building MFEM](https://mfem.org/building/)

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
   sudo make install  # optional
   ```

3. Build and install HYPRE

   ```bash
   cd <mypath>
   cd hypre/src && \
   mkdir build && \
   cd build && \
   cmake .. && \
   make -j $(nproc --all) && \
   sudo make install # optional
   ```

4. Build and install MFEM

   ```bash
   cd <mypath>
   cd mfem && \
   mkdir build && \
   cd build && \
   cmake .. && \
   make -j $(nproc --all) && \
   sudo make install # optional
   ```
