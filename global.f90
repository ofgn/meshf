! ------------------------------------------------------------------------------
! @file global.f90
! @brief Global module for shared variables and functions.
! @date 2024-07-01
! ------------------------------------------------------------------------------
module global
    use iso_fortran_env, only: int16, int32, int64, real64, real128, &
        logical_kinds, output_unit, error_unit

    implicit none

    integer(int32), parameter :: BOOL8 = logical_kinds(1)                       ! 8-bit boolean

    integer(int32), parameter :: PARALLEL_THRESHOLD = 10**5                     ! When to switch to parallel code

    integer(int16), parameter :: STDOUT = output_unit                           ! Standard output
    integer(int16), parameter :: STDERR = error_unit                            ! Standard error

    integer(int16), parameter :: LEVEL_DEBUG = 0                                ! Debug level
    integer(int16), parameter :: LEVEL_INFO = 1                                 ! Info level
    integer(int16), parameter :: LEVEL_WARNING = 2                              ! Warning level
    integer(int16), parameter :: LEVEL_ERROR = 3                                ! Error level

    integer(int16) :: log_unit = STDOUT                                         ! Default log unit
    logical(logical_kinds(1)) :: log_enabled = .true.                           ! Enable or disable logging
    integer(int16) :: log_level = LEVEL_INFO                                    ! Default log level

contains

    ! --------------------------------------------------------------------------
    ! @brief Configure the logging system.
    ! @param[in] unit File unit to use for logging (optional).
    ! @param[in] enabled Logical flag to enable/disable logging (optional).
    ! @param[in] level Minimum log level to display (optional).
    ! --------------------------------------------------------------------------
    subroutine configure_logging(unit, enabled, level)
        implicit none

        integer(int16), intent(in), optional :: unit                            ! File unit for logging
        logical(logical_kinds(1)), intent(in), optional :: enabled              ! Enable/disable logging
        integer(int16), intent(in), optional :: level                                  ! Minimum log level

        if (present(unit)) then
            log_unit = unit
        end if
        if (present(enabled)) then
            log_enabled = enabled
        end if
        if (present(level)) then
            log_level = level
        end if
    end subroutine configure_logging

    ! --------------------------------------------------------------------------
    ! @brief Reset logging to default settings.
    ! --------------------------------------------------------------------------
    subroutine reset_logging()
        implicit none

        log_unit = STDOUT                                                       ! Reset to standard output
        log_enabled = .true.                                                    ! Enable logging by default
        log_level = LEVEL_INFO                                                  ! Default to info level
    end subroutine reset_logging

    ! --------------------------------------------------------------------------
    ! @brief Log a message at a specified level.
    ! @param[in] message Message to log.
    ! @param[in] level Level of the message (0 = debug, 1 = info, etc.).
    ! --------------------------------------------------------------------------
    subroutine report(message, level)
        implicit none

        character(len=*), intent(in) :: message                                 ! Message to log
        integer(int16), intent(in) :: level                                     ! Log level of the message

        character(len=10) :: prefix                                             ! Prefix for log message

        if (.not. log_enabled .or. level < log_level) then
            return
        end if

        select case (level)
            case (LEVEL_DEBUG)
                prefix = "[DEBUG]   "
            case (LEVEL_INFO)
                prefix = "[INFO]    "
            case (LEVEL_WARNING)
                prefix = "[WARNING] "
            case (LEVEL_ERROR)
                prefix = "[ERROR]   "
            case default
                prefix = "[UNKNOWN] "
        end select

        write(log_unit, '(A)') prefix // trim(message)
    end subroutine report
end module global
