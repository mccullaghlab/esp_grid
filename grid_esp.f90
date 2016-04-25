
!requires lapack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Modules !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module inputData

        integer deltaStep
        integer nThreads

endmodule inputData

module atomData

        integer nAtoms
        real (kind=8), allocatable :: atomPos(:,:)
        real (kind=8), allocatable :: atomCharges(:)
        character*80 atomPsfFile
        character*80, allocatable :: atomDcdFiles(:)
        integer nDcdFiles
        character*80 outFile
        real (kind=8), allocatable :: atomEspMat(:)

endmodule atomData

module cgData

        integer nCgs
        real (kind=8), allocatable :: cgPos(:,:)
        real (kind=8), allocatable :: cgCharges(:)
        character*80, allocatable :: cgDcdFiles(:)
        real (kind=8), allocatable :: cgEspMat(:,:)

endmodule cgData

module gridData

        integer nGrids(3)
        integer totalGrids
        real (kind=8) gridMin(3)
        real (kind=8) gridMax(3)
        real (kind=8) deltaGrid
        real (kind=8) cutoff2
        integer gridCut
!        real (kind=8), save :: maxThresh2 = 400.0
        real (kind=8), save :: minThresh2 = 0.01

endmodule gridData

module integralData

        real (kind=8), allocatable :: A(:,:)
        real (kind=8), allocatable :: B(:,:)
        real (kind=8), allocatable :: C(:,:)
        real (kind=8), allocatable :: intCgCharges(:)
        real (kind=8) gridRss
        real (kind=8) intRss

endmodule integralData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Main Program !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program esp_grid
        use atomData
        use cgData
        use inputData
        use gridData, only : totalGrids
        implicit none
        character*80 cfgFile
        real (kind=8) omp_get_wtime
        real (kind=8) ti,tf

        ti = omp_get_wtime()
        ! read config file, number of openMP threads and stride from command line
        call parse_command_line(cfgFile,nThreads,deltaStep)
        
        ! read the rest of the config parameters from the config file
        call parse_config_file(cfgFile)

        ! read the atom psf file to obtain number of atoms and charges.  atomPos and charge array are allocated in this routine
        call read_psf_file()

        ! allocate some grid and integral arrays
        call allocate_arrays()

        ! read the trajectories and perform the fits
        call read_trajectory_compute_grid_esp()

        tf = omp_get_wtime()
        write(*,'("Total time elapsed:",f8.3)') tf-ti

endprogram esp_grid



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Subroutines  !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! allocate arrays
subroutine allocate_arrays()
        use atomData
        use gridData
        use cgData
        use integralData
        implicit none

        ! grid based array containing the ESP multiplied by the atomic charges
        allocate(atomEspMat(totalGrids))
        ! grid based array containing 1/dist between grids and CG sites
        allocate(cgEspMat(totalGrids,nCgs))
        ! zero out each array as they will be added to
        atomEspMat = 0.0
        cgEspMat = 0.0

        ! integral based matrices
        allocate(A(nCgs,nCgs),B(nCgs,nAtoms),C(nAtoms,nAtoms),intCgCharges(nCgs))

        A = 0.0
        B = 0.0
        C = 0.0

endsubroutine allocate_arrays

subroutine read_trajectory_compute_grid_esp()
        use inputData
        use atomData
        use cgData
        use gridData
        use integralData
        implicit none
        integer atom1
        integer k
        integer nSteps
        integer step
        real (kind=8) centerVec(3)
        real (kind=8) gridCenter(3)
        real (kind=8) temp
        integer filenum
        integer totalSteps

        gridCenter = (gridMax+gridMin) / 2.0
        print*, "Grid center:", gridCenter(1),gridCenter(2),gridCenter(3)
        totalSteps = 0
        do filenum = 1, nDcdFiles
                call read_dcd_header(atomDcdFiles(filenum),nAtoms,nSteps,20)
                call read_dcd_header(cgDcdFiles(filenum),nCgs,nSteps,30)
                write(*,'("Starting step loop for ",i10," number of steps")') nSteps
                
                do step=1,nSteps
                        if (mod(step,deltaStep)==0) then
                                totalSteps = totalSteps + 1
                                write(*,'("Reading step ",i10," of ",i10)') step, nSteps
                                call read_dcd_step(atomPos,nAtoms,20)
                                call read_dcd_step(cgPos,nCgs,30)
                                write(*,'("read trajectory data")') 
                                ! make sure molecule is in center of grid?
                                do k=1,3
                                        centerVec(k) = gridCenter(k) - sum(atomPos(:,k)) / dble(nAtoms)
                                        atomPos(:,k) = atomPos(:,k) + centerVec(k)
                                        cgPos(:,k) = cgPos(:,k) + centerVec(k)
                                enddo     
                                print*, "sum of atom x component:", sum(atomPos(:,1))     
                                print*, "sum of cg x component:", sum(cgPos(:,1))     
                            
                                ! compute ESP matrices.  Atomic one is distance to grid points times atomic charge.  
                                ! CG is only distance to grid points.
                                call compute_atom_grid_esp(atomPos,nAtoms,atomCharges,atomEspMat)
                                call compute_cg_grid_esp(cgPos,nCgs,cgEspMat)
                                ! fit - do we need to fit every time and average the CG charges at the end

                                ! update A, B and C matrices for integral method
                                call update_A_B_C_matrices(atomPos,nAtoms,cgPos,nCgs,A,B,C)
                        endif
                enddo
                close(20)
                close(30)
        enddo
        print*, "done reading data!"
        !average
        atomEspMat = atomEspMat / dble(totalSteps)
        cgEspMat = cgEspMat / dble(totalSteps)

        print*, "fitting grid charges"
        ! fit?
        call fit_cg_charges(atomEspMat,nAtoms,cgEspMat,nCgs,totalGrids,cgCharges)

        print*, "fitting integral charges"
        ! fit using integral approach
        call integral_fit_charges(A, B, C, atomCharges, intCgCharges, nAtoms, nCgs, intRss)

        ! print CG charges
        open(35,file=outFile)
        write(35,'(a10,a20,a20,a20)') "CG Site", "Grid Charges", "Integral Charges", "Squared difference"
        write(*,'(a10,a20,a20,a20)') "CG Site", "Grid Charges", "Integral Charges", "Squared difference"
        temp = 0
        do k=1, nCgs
                temp = temp + (cgCharges(k)-intCgCharges(k))**2
                write(35,'(i10,f20.10,f20.10,f20.10)') k , cgCharges(k), intCgCharges(k), (cgCharges(k)-intCgCharges(k))**2
                write(*,'(i10,f20.10,f20.10,f20.10)') k , cgCharges(k), intCgCharges(k), (cgCharges(k)-intCgCharges(k))**2
        enddo    
        write(35,'("Sums:", a20,a20,a20,a20)') "Atomic", "Grid Charges", "Integral Charges", "Sum of sq diff" 
        write(35,'(4x,f20.10,f20.10,f20.10,f20.10)') sum(atomCharges), sum(cgCharges), sum(intCgCharges), temp
        close(35)
        write(*,'("Sums:", a20,a20,a20,a20)') "Atomic", "Grid Charges", "Integral Charges", "Sum of sq diff"
        write(*,'(4x,f20.10,f20.10,f20.10,f20.10)') sum(atomCharges), sum(cgCharges), sum(intCgCharges), temp
        close(35)
!        write(*,'("CG grid charge:", f10.5, "CG int charges:", f10.5, " total atom charge:", f10.5)') sum(cgCharges), sum(intCgCharges), sum(atomCharges)

endsubroutine read_trajectory_compute_grid_esp

! solve the least squares problem to fit CG charges to all-atom ESP 
subroutine fit_cg_charges(atomEspMat,nAtoms,cgEspMat,nCgs,totalGrids,cgCharges)
        implicit none
        integer nCgs
        integer nAtoms
        integer totalGrids
        real (kind=8) atomEspMat(totalGrids,1)
        real (kind=8) cgEspMat(totalGrids,nCgs)
        real (kind=8) cgCharges(nCgs)
        integer info
        integer lwork 
        real (kind=8) workQuery(1)
        real (kind=8), allocatable :: work(:)

        print*, "fitting charges..."
        lwork = -1
        call dgels('N', totalGrids, nCgs, 1, cgEspMat, totalGrids, atomEspMat, totalGrids, workQuery, lwork, info)
        ! formula for optimal lwork value
!        lwork = max(1, min(totalGrids,nCgs) + max(min(totalGrids,nCgs),1)*8)
        lwork = int(workQuery(1))
        print*, "lwork = ", lwork
        allocate(work(lwork))

        call dgels('N', totalGrids, nCgs, 1, cgEspMat, totalGrids, atomEspMat, totalGrids, work, lwork, info)

        print*, "solved linear least squares problem"

        cgCharges = atomEspMat(1:nCgs,1)

endsubroutine fit_cg_charges

!
subroutine compute_cg_grid_esp(cgPos,nCgs,gridEspMat)
        use inputData, only : nThreads
        use gridData
        implicit none
        integer nCgs
        real (kind=8) cgPos(nCgs,3)
        real (kind=8) gridEspMat(totalGrids,nCgs)
        integer atom
        real (kind=8) diff(3)
        real (kind=8) dist2, dist
        integer x, y, z, gridCount(2), gridCountSum, k
        integer gridStart(3), gridStop(3), gridAtom(3)

        !$omp parallel num_threads(nThreads) private(atom,gridStart,x,y,z,gridCount,diff,dist2,dist) shared(cutoff2,cgPos,nCgs,nGrids,gridEspMat,gridMin,deltaGrid,minThresh2)
        !$omp do
        do atom = 1, nCgs
                gridAtom = int( (cgPos(atom,:)-gridMin) / deltaGrid ) + 1
                ! determine grid cutoffs based on distance cutoff
                do k=1,3
                        if ( gridAtom(k) > gridCut) then
                                gridStart(k) = gridAtom(k) - gridCut
                        else
                                gridStart(k) = 1
                        endif
                        if ( gridAtom(k) < nGrids(k) - gridCut) then
                                gridStop(k) = gridAtom(k) + gridCut
                        else
                                gridStop(k) = nGrids(k)
                        endif
                enddo
                ! loop through all possible grid points
                do x=gridStart(1),gridStop(1)
                        diff(1) = cgPos(atom,1) - ( gridMin(1) + (x-1)*deltaGrid )
                        gridCount(1) = (x-1) * nGrids(2)*nGrids(3) 
                        do y = gridStart(2), gridStop(2)
                                diff(2) = cgPos(atom,2) - ( gridMin(2) + (y-1)*deltaGrid ) 
                                gridCount(2) = (y-1) * nGrids(3)
                                do z= gridStart(3), gridStop(3)
                                        diff(3) = cgPos(atom,3) - ( gridMin(3) + (z-1)*deltaGrid ) 
                                        gridCountSum = sum(gridCount) + z
                                        dist2 = dot_product(diff,diff)
                                        ! check to see if within cutoffs
                                        if (dist2 > minThresh2 .and. dist2 < cutoff2) then
                                                dist = sqrt(dist2)
                                                gridEspMat(gridCountSum,atom) = gridEspMat(gridCountSum,atom) + 1.0/dist
                                        endif
                                enddo

                        enddo
                enddo
        enddo
        !$omp end do
        !$omp end parallel

endsubroutine compute_cg_grid_esp

!compute the center of geometry positions of the BPs
subroutine compute_atom_grid_esp(atomPos,nAtoms,atomCharges,gridEspMat)
        use inputData, only : nThreads
        use gridData
        implicit none
        integer nAtoms
        real (kind=8) atomPos(nAtoms,3)
        real (kind=8) atomCharges(nAtoms)
        real (kind=8) gridEspMat(totalGrids)
        integer atom
        real (kind=8) diff(3)
        real (kind=8) dist2, dist
        integer x, y, z, gridCount(2), gridCountSum, k
        integer gridStart(3), gridStop(3), gridAtom(3)

        !$omp parallel num_threads(nThreads) private(atom,gridStart,x,y,z,gridCount,diff,dist2,dist) shared(cutoff2,atomPos,nAtoms,nGrids,gridEspMat,gridMin,deltaGrid,minThresh2)
        !$omp do
        do atom = 1, nAtoms
!                write(*,'("atom ", i5, " of ", i5)') atom, nAtoms
                ! determine grid point for current atom position
                gridAtom = int( (atomPos(atom,:)-gridMin) / deltaGrid ) + 1
                ! determine grid cutoffs based on distance cutoff
                do k=1,3
                        if ( gridAtom(k) > gridCut) then
                                gridStart(k) = gridAtom(k) - gridCut
                        else
                                gridStart(k) = 1
                        endif
                        if ( gridAtom(k) < nGrids(k) - gridCut) then
                                gridStop(k) = gridAtom(k) + gridCut
                        else
                                gridStop(k) = nGrids(k)
                        endif
                enddo
                ! loop through all possible grid points
                do x=gridStart(1),gridStop(1)
                        diff(1) = atomPos(atom,1) - ( gridMin(1) + (x-1)*deltaGrid )
                        gridCount(1) = (x-1) * nGrids(2)*nGrids(3) 
                        do y = gridStart(2), gridStop(2)
                                diff(2) = atomPos(atom,2) - ( gridMin(2) + (y-1)*deltaGrid ) 
                                gridCount(2) = (y-1) * nGrids(3)
                                do z= gridStart(3), gridStop(3)
                                        diff(3) = atomPos(atom,3) - ( gridMin(3) + (z-1)*deltaGrid ) 
                                        gridCountSum = sum(gridCount) + z
                                        dist2 = dot_product(diff,diff)
                                        if (dist2 > minThresh2 .and. dist2 < cutoff2) then
                                                dist = sqrt(dist2)
                                                gridEspMat(gridCountSum) = gridEspMat(gridCountSum) + atomCharges(atom)/dist
                                        endif
                                enddo

                        enddo
                enddo
        enddo
        !$omp end do
        !$omp end parallel

endsubroutine compute_atom_grid_esp

subroutine parse_command_line(cfgFile,nThreads,deltaStep)
        implicit none
        integer i
        character*30 arg
        character*80 cfgFile
        integer nThreads
        integer deltaStep
        logical deltaStepFlag
        logical cfgFileFlag
        logical nThreadsFlag

        deltaStepFlag = .false.
        nThreadsFlag = .false.
        cfgFileFlag = .false.

        i=1
        do 
   
                call get_command_argument(i, arg) 
   
                select case (arg) 
   
                case ('-cfg')
                        i = i+1
                        call get_command_argument(i,cfgFile)
                        cfgFileFlag=.true.
                        print*, "CFG File: ", cfgFile
                case ('-stride') 
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(i10)') deltaStep
                        print*, "Stride: ", deltaStep
                        deltaStepFlag = .true.
                case ('-np') 
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(i10)') nThreads
                        print*, "Number of threads: ", nThreads
                        nThreadsFlag = .true.
                case default 
                        print '(a,a,/)', 'Unrecognized command-line option: ', arg 
                        print*, 'Usage: esp_grid.x -cfg [cfg file] -stride [delta step size] -np [number of threads]'
                        stop 
                end select 
                i = i+1
                if (i.ge.command_argument_count()) exit
        enddo

        if (cfgFileFlag.eqv..false.) then
                write(*,'("Must provide a cfg file using command line argument -cfg [cfg file name]")')
                stop
        endif
        if (deltaStepFlag.eqv..false.) then
                write(*,'("Using default step size of 1.  Change this with command line argument -stride [delta step size]")')
                deltaStep = 1
        endif
        if (nThreadsFlag.eqv..false.) then
                write(*,'("Using default of 1 thread.  Change this with command line argument -np [number of threads]")')
                nThreads = 1
        endif

endsubroutine parse_command_line


!Read atomic charges from some file
subroutine read_psf_file
        use atomData
        implicit none
!        real (kind=8), parameter :: chargeConvert = 553.43949573 ! to convert to kT/e with dielectric of 1
        real (kind=8), parameter :: chargeConvert = 1.0 ! 
        integer atom, j
        character*6 check           !character to check if NATOM is in the line
        character*8 numChar         !character to read number of atoms.  must be converted to an integer
        character*4 atomCheck
        character*24 posChar
        integer ios

        !open the psf file
        print*, "psf file name here:", atomPsfFile
        open(10,file=atomPsfFile)

        !run through the header of the file and look for the number of atoms
        do 

                read(10,'(a8,2x,a6)') numChar, check

                !if we read the number of atoms exit this do loop
                if (check.eq.'NATOM ') then
                        !this converts the character natoms_char to the integer natoms
                        read(numChar,*) nAtoms
                        write(*,*) "Number of atoms=", nAtoms
                        !Now that we know the number of atoms we must allocate the arrays
                        allocate(atomCharges(nAtoms))
                        allocate(atomPos(nAtoms,3))
                        !Now we loop through the number of atoms and read the pertinent information
                        do atom=1,nAtoms

                                read(10,'(34x,f10.6)') atomCharges(atom) 
                                atomCharges(atom) = atomCharges(atom) * chargeConvert
                                !add to total charge

                        enddo
                        write(*,*) "Total charge in atom psf=", sum(atomCharges)
                elseif (check.eq.'NBOND:') then
                        exit
                endif

        enddo

        close(10)

endsubroutine read_psf_file


!read config information from file
subroutine parse_config_file(cfgFile)
        use gridData
        use atomData
        use cgData
        implicit none
        character*80 cfgFile
        character*200 line
        character*80,firstWord
        integer i, k
        character*30 arg
        character*30 sep
        real (kind=8) cutoff
        integer ios
        integer atomDcdFileCount
        integer cgDcdFileCount
        logical atomPsfFlag
        logical atomDcdFlag
        logical cgDcdFlag
        logical outFileFlag
        logical gridMaxFlag
        logical gridMinFlag
        logical deltaGridFlag
        logical nCgFlag
        logical cutoffFlag

        atomPsfFlag = .false.
        atomDcdFlag = .false.
        cgDcdFlag = .false.
        outFileFlag = .false.
        gridMaxFlag = .false.
        gridMinFlag = .false.
        deltaGridFlag = .false.
        nCgFlag = .false.
        cutoffFlag = .false.

        open(12,file=cfgFile)
        atomDcdFileCount = 1
        cgDcdFileCount = 1
        do while(ios>=0)
                read(12,'(a600)',IOSTAT=ios) line
                call split(line,'=',firstWord, sep)
                if (line .ne. "") then
                        if (firstWord .eq. "numdcd") then
                                read(line,'(i10)') nDcdFiles
                                write(*,*) "Number of dcd files: ", nDcdFiles
                                allocate(atomDcdFiles(nDcdFiles),cgDcdFiles(nDcdFiles))
                                atomDcdFlag = .true.
                                cgDcdFlag = .true.
                        else if (firstWord .eq. "atomdcdfile" .and. atomDcdFlag .eqv. .true.) then
                                atomDcdFiles(atomDcdFileCount) = line
                                write(*,*) "Dcd file number ", atomDcdFileCount,":", atomDcdFiles(atomDcdFileCount)
                                atomDcdFileCount = atomDcdFileCount + 1
                        else if (firstWord .eq. "cgdcdfile" .and. cgDcdFlag .eqv. .true.) then
                                cgDcdFiles(cgDcdFileCount) = line
                                write(*,*) "Dcd file number ", cgDcdFileCount,":", cgDcdFiles(cgDcdFileCount)
                                cgDcdFileCount = cgDcdFileCount + 1
                        else if (firstWord .eq. "psffile") then
                                atomPsfFile = line
                                write(*,*) "psf file:", atomPsfFile
                                atomPsfFlag = .true.
                        else if (firstWord .eq. "outfile") then
                                outFile = line
                                write(*,*) "out file:", outFile
                                outFileFlag = .true.
                        else if (firstWord .eq. "deltagrid") then
                                read(line,'(f10.5)') deltaGrid
                                write(*,*) "deltaGrid:", deltaGrid
                                deltaGridFlag = .true.
                        else if (firstWord .eq. "cutoff") then
                                read(line,'(f10.5)') cutoff
                                cutoff2 = cutoff * cutoff
                                write(*,*) "cutoff:", cutoff
                                cutoffFlag = .true.
                        else if (firstWord .eq. "ncg") then
                                read(line,'(i10)') nCgs
                                write(*,*) "nCgs:", nCgs
                                allocate(cgPos(nCgs,3))
                                nCgFlag = .true.
                        else if (firstWord .eq. "maxgrid") then
                                do k = 1,3
                                        call split(line,' ',firstWord, sep)
                                        read(firstWord,'(f10.5)') gridMax(k)
                                enddo
                                write(*,*) "max grid:", gridMax(1), gridMax(2), gridMax(3)
                                gridMaxFlag = .true.
                        else if (firstWord .eq. "mingrid") then
                                do k = 1,3
                                        call split(line,' ',firstWord, sep)
                                        read(firstWord,'(f10.5)') gridMin(k)
                                enddo
                                write(*,*) "min grid:", gridMin(1),gridMin(2),gridMin(3)
                                gridMinFlag = .true.
                        endif
                endif
        enddo

        close(12)


        if (atomPsfFlag.eqv..false.) then
                write(*,'("Must provide a psf file using command line argument -psf [psf file name]")')
                stop
        endif
        if (nCgFlag.eqv..false.) then
                write(*,'("Must provide the number of cg sites.  Use ncg = [number of cg sites] in config file")')
                stop
        endif
        if (cutoffFlag.eqv..false.) then
                write(*,'("cutoff not defined.  using default of 20.0.  change with cutoff = [cutoff] in config file")')
                cutoff = 20.0
                cutoff2 = cutoff * cutoff
        endif
        if (atomDcdFlag.eqv..false.) then
                write(*,'("Must provide a atom dcd file using command line argument -adcd [atom dcd file name]")')
                stop
        endif
        if (outFileFlag.eqv..false.) then
                write(*,'("Must provide an output file name  using command line argument -o [output file name]")')
                stop
        endif
        if (gridMaxFlag.eqv..false.) then
                write(*,'("Using default max radius of 400. Change this with command line argument -max [max radius]")')
                gridMax = 100.0
        endif
        if (gridMinFlag.eqv..false.) then
                write(*,'("Using default min radius of 0. Change this with command line argument -min [min radius]")')
                gridMin = 0.0
        endif
        if (deltaGridFlag.eqv..false.) then
                write(*,'("Using default of deltaGrid=1.0.")')
                deltaGrid=1.0
        endif
        totalGrids = 1.0
        do k=1,3
                nGrids(k) = int( (gridMax(k)-gridMin(k))/deltaGrid) + 1
                totalGrids = totalGrids*nGrids(k)
        enddo
        allocate(cgCharges(nCgs))
        cgCharges = 1.0
        gridCut = int(cutoff/deltaGrid) + 1

endsubroutine parse_config_file


!subroutine to compute CG charges from input matrices A and B.  The residual sum of squares is calculated with aid of all-atom matrix C
subroutine integral_fit_charges(A, B, C, atomCharges, cgCharges, nAtoms, nCg, rss)
        implicit none
        integer nAtoms
        integer nCg
        real (kind=8) A(nCg,nCg)
        real (kind=8) ATemp(nCg,nCg)
        real (kind=8) D(nCg-1,nCg)
        real (kind=8) B(nCg,nAtoms)
        real (kind=8) C(nAtoms,nAtoms)
        real (kind=8) BTemp(nCg,nAtoms)
        real (kind=8) cgCharges(nCg,1)
        real (kind=8) atomCharges(nAtoms)
        real (kind=8) atomChargesM(nAtoms,1)
        real (kind=8) newB(nCg)
        real (kind=8) temp(1,1)
        real (kind=8) rss
        integer j, i, k
        !lapack routine variables
        integer ipiv(nCg)
        integer info

        !First we need to modify A and B to have the correct matrix properties 
        !create D matrices
        D=0.0
        do i=1,nCg-1
                D(i,i)=1.0
                D(i,i+1)=-1.0
        enddo 
        !multiply A by D0 and B by D1 giving new matrices A and B the correct behavior
        ATemp(1:(nCg-1),:) = matmul(D,A)
        BTemp(1:(nCg-1),:) = matmul(D,B)
        !generate new matrices with last line having 1.0s forcing charge conservation
        ATemp(nCg,:) = 1.0
        BTemp(nCg,:) = 1.0

        ! multiple right hand side of equation by atomic charges
        newB = matmul(BTemp,atomCharges)
        ! determine the solution to system of linear equations ATemp*X = newB
        call dgesv(nCg,1, ATemp,nCg,ipiv,newB,nCg,info)
        
        cgCharges(:,1) = real(newB(1:nCg))
        atomChargesM(:,1) = atomCharges

        ! compute residual sum of squares
        temp = matmul(transpose(atomChargesM),matmul(C,atomChargesM))+matmul(transpose(cgCharges),matmul(A,cgCharges))-2*matmul(transpose(cgCharges),matmul(B,atomChargesM))
        rss = temp(1,1)

endsubroutine integral_fit_charges

!requires lapack
subroutine update_A_B_C_matrices(atomPos,nAtoms,cgPos,nCgs,A,B,C)
        implicit none
        integer nAtoms
        integer nCgs
        real (kind=8) atomPos(nAtoms,3)
        real (kind=8) atomCharges(nAtoms)
        real (kind=8) cgPos(nCgs,3)
        real (kind=8) A(nCgs,nCgs)
        real (kind=8) B(nCgs,nAtoms)
        real (kind=8) C(nAtoms,nAtoms)
        real (kind=8) dist, temp
        !loop indeces
        integer cgSite1, cgSite2
        integer j
        integer atom1, atom2

        ! populate A matrix with negative distances between CG sites
        do cgSite1 = 1, nCgs-1
                do cgSite2 = cgSite1+1,nCgs
                        dist = 0
                        do j=1,3
                                temp = cgPos(cgSite1,j)-cgPos(cgSite2,j)
                                dist = dist + temp*temp
                        enddo
                        dist = sqrt(dist)
                        A(cgSite1,cgSite2) = A(cgSite1,cgSite2)-dist
                        !symmetrize the matrix
                        A(cgSite2,cgSite1) = A(cgSite1,cgSite2)
                enddo
        enddo

        ! populate B matrix with negative distance between CG sites and atoms
        do cgSite1 = 1, nCgs
                do atom1 = 1,nAtoms
                        dist = 0
                        do j=1,3
                                temp = cgPos(cgSite1,j)-atomPos(atom1,j)
                                dist = dist + temp*temp
                        enddo
                        B(cgSite1,atom1) = B(cgSite1,atom1)-sqrt(dist)
                enddo
        enddo

        ! populate C matrix with negative distance between atoms
        do atom1 = 1, nAtoms-1
                do atom2 = atom1+1, nAtoms
                        dist = 0
                        do j=1,3
                                temp = atomPos(atom1,j)-atomPos(atom2,j)
                                dist = dist + temp*temp
                        enddo
                        C(atom1,atom2) = C(atom1,atom2)-sqrt(dist)
                        ! symmetrize the matrix
                        C(atom2,atom1) = C(atom1,atom2)
                enddo
        enddo

endsubroutine update_A_B_C_matrices


