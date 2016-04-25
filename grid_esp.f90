
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
        call parse_command_line(cfgFile,nThreads,deltaStep)
        
        call parse_config_file(cfgFile)

        call read_psf_file()

        allocate(atomEspMat(totalGrids))
        allocate(cgEspMat(totalGrids,nCgs))
        atomEspMat = 0.0
        cgEspMat = 0.0

        call read_trajectory_compute_grid_esp()

        tf = omp_get_wtime()
        write(*,'("Total time elapsed:",f8.3)') tf-ti

endprogram esp_grid



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Subroutines  !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_trajectory_compute_grid_esp()
        use inputData
        use atomData
        use cgData
        use gridData
        implicit none
        real (kind=8), parameter :: pi = 3.1415926535
        integer atom1
        integer k
        integer nSteps
        integer step
        real (kind=8) centerVec(3)
        real (kind=8) gridCenter(3)
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
                            
                                ! compute the center of geometry positions of the BPs
                                call compute_atom_grid_esp(atomPos,nAtoms,atomCharges,atomEspMat)
                                call compute_cg_grid_esp(cgPos,nCgs,cgEspMat)
                                ! fit?

                        endif
                enddo
                close(20)
                close(30)
        enddo
        print*, "done reading data!"
        !average
        atomEspMat = atomEspMat / dble(totalSteps)
        cgEspMat = cgEspMat / dble(totalSteps)

        ! fit?
        call fit_cg_charges(atomEspMat,nAtoms,cgEspMat,nCgs,totalGrids,cgCharges)

        ! print CG charges
        open(35,file=outFile)
        do k=1, nCgs
                write(35,'(i10,f20.10)') k , cgCharges(k)
                write(*,'(i10,f20.10)') k , cgCharges(k)
        enddo        
        close(35)
        write(*,'("total CG charge:", f10.5, " total atom charge:", f10.5)') sum(cgCharges), sum(atomCharges)

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


