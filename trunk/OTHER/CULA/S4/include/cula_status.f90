
! types
module cula_type

    implicit none

    ! culastatus
    integer, parameter :: culanoerror = 0
    integer, parameter :: culanotinitialized = 1
    integer, parameter :: culanohardware = 2
    integer, parameter :: culainsufficientruntime = 3
    integer, parameter :: culainsufficientcomputecapability = 4
    integer, parameter :: culainsufficientmemory = 5
    integer, parameter :: culafeaturenotimplemented = 6
    integer, parameter :: culaargumenterror = 7
    integer, parameter :: culadataerror = 8
    integer, parameter :: culablaserror = 9
    integer, parameter :: cularuntimeerror = 10
    integer, parameter :: culabadstorageformat = 11
    integer, parameter :: culaunspecifiederror = 12

end module cula_type

module cula_status

    use cula_type

    implicit none

    interface
        integer function cula_initialize()
        end function
    end interface

    interface
        subroutine cula_shutdown()
        end subroutine
    end interface

    interface
      integer function cula_get_error_info_string(status,info,buf,bufsize)
        integer, intent(in) :: status
        integer, intent(in) :: info
        character*(*), intent(out) :: buf
        integer, intent(in) :: bufsize
      end function
    end interface

    interface
        integer function cula_get_last_status()
        end function
    end interface

    interface
        integer function cula_get_error_info()
        end function
    end interface

    interface
        subroutine cula_free_buffers()
        end subroutine
    end interface

    interface
        integer function cula_get_version()
        end function
    end interface

    interface
        integer function cula_get_cuda_minimum_version()
        end function
    end interface

    interface
        integer function cula_get_cuda_runtime_version()
        end function
    end interface

    interface
        integer function cula_get_cuda_driver_version()
        end function
    end interface

    interface
        integer function cula_get_cublas_minimum_version()
        end function
    end interface

    interface
        integer function cula_get_cublas_runtime_version()
        end function
    end interface

contains

    subroutine cula_check_status(status)
        integer :: status
        integer :: s
        integer :: info
        integer :: culaversion
        integer, parameter :: bufsize = 256
        character*256 :: buf

        if (status .ne. culanoerror) then
            culaversion = cula_get_version()
            info = cula_get_error_info()

            ! There is an error in CULA R14's Fortran cula_get_error_info_string
            ! For version R14 and earlier
            ! So we have implemented a lightweight version here
            if ( culaversion .le. 14000 ) then
                if (status .eq. 1) then
                    write(*,*) 'CULA is not initialized'
                else if (status .eq. 2) then
                    write(*,*) 'no hardware has been detected'
                else if (status .eq. 3) then
                    write(*,*) 'insufficient CUDA runtime was loaded'
                else if (status .eq. 4) then
                    write(*,*) 'hardware lacks sufficient compute capability for requested operation'
                else if (status .eq. 5) then
                    write(*,*) 'no hardware has been detected'
                else if (status .eq. 6) then
                    write(*,*) 'the CULA function invoked is not implemented'
                else if (status .eq. 7) then
                    write(*,*) 'invalid value for parameter ', info
                else if (status .eq. 8) then
                    write(*,*) 'data error (', info ,')'
                else if (status .eq. 9) then
                    write(*,*) 'blas error (', info ,')'
                else if (status .eq. 10) then
                    write(*,*) 'runtime error - consult the CUDA error list for code ', info
                else
                    write(*,*) 'an unknown error has occurred'
                endif
            else
                s = cula_get_error_info_string(status, info, buf, bufsize)
                write(*,'(a)') buf
            endif
            if (status .eq. culainsufficientcomputecapability) then
                stop 0
            else
                stop 1
            endif
        endif
    end subroutine

end module cula_status

