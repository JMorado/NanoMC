MODULE monte_carlo_module
    USE simulation_module
    USE energy_module
    USE pbc_module
    USE displacement_module

    IMPLICIT NONE
    ABSTRACT INTERFACE
        FUNCTION energy_function(simulation_instance, index, coord)
            IMPORT                            :: Simulation
            REAL(8)                           :: energy_function
            TYPE(Simulation), INTENT(IN)      :: simulation_instance
            REAL(8),           INTENT(IN)     :: coord(:)
            INTEGER(4),        INTENT(IN)     :: index
        END FUNCTION energy_function
    END INTERFACE


CONTAINS
    SUBROUTINE mc_move(simulation_instance, energy_func, accepted, o)
        !========================================================================!
        ! This subroutine attempts to perform a standard MC move by displacing   !
        ! a particle.                                                            !
        !                                                                        !
        ! Reference:                                                             !
        ! Understanding Molecular Simulation                                     !
        ! Smit B and Frenkel D                                                   !
        ! Elsevier, 2002. https://doi.org/10.1016/B978-0-12-267351-1.X5000-7.    !
        ! pp. 32-33                                                              !
        !------------------------------------------------------------------------!
        ! simulation_instance      (inout) : simulation instance                 !
        ! energy_func                 (in) : pointer to energy function          !
        ! accepted                    (in) : flag that signal move acceptance    !
        ! o                           (in) : particle index                      !
        ! displacement_type           (in) : "gaussian" or "uniform"             !
        !========================================================================!
        IMPLICIT NONE
        TYPE(Simulation), INTENT(INOUT)     :: simulation_instance
        LOGICAL(1), INTENT(OUT)             :: accepted
        INTEGER(4), INTENT(IN)              :: o
        REAL(8)                             :: deltaE, e_old, e_new, prob_acceptance
        REAL(8)                             :: coord_tmp(3)
        REAL(8)                             :: dr(3)
        PROCEDURE(energy_function), POINTER :: energy_func

        accepted = .FALSE.

        ! Energy old configuration
        e_old = energy_func(simulation_instance, o, simulation_instance%coord(:,o))

        ! Generate new configuration
        ! Choose correct energy functions
        SELECT CASE (simulation_instance%displacement_type)
        CASE("gaussian")
            dr =  get_displacement_gaussian( simulation_instance%disp )
        CASE("uniform")
            dr =  get_displacement_uniform( simulation_instance%disp )
        END SELECT

        coord_tmp = simulation_instance%coord(:,o) + dr(:)

        ! Apply PBC
        CALL pbc(coord_tmp, (/simulation_instance%cell_instance%a,simulation_instance%cell_instance%b,&
                simulation_instance%cell_instance%c/))

        ! Energy new configuration
        e_new = energy_func(simulation_instance, o, coord_tmp(:))

        ! Energy difference
        deltaE = e_new-e_old

        ! Metropolis criteria
        IF (deltaE .le. 0) THEN
            simulation_instance%dr = dr
            simulation_instance%deltaE = deltaE

            simulation_instance%coord(:,o) = coord_tmp
            simulation_instance%coord_eff(:,o) =  simulation_instance%coord_eff(:,o) + dr

            accepted = .TRUE.
        ELSE IF ((deltaE / simulation_instance%temperature) .le. 13) THEN
            prob_acceptance = min(exp(- deltaE / simulation_instance%temperature),1.0)

            IF (prob_acceptance .ge. RAND()) THEN
                simulation_instance%dr = dr
                simulation_instance%deltaE = deltaE

                simulation_instance%coord(:,o) = coord_tmp
                simulation_instance%coord_eff(:,o) =  simulation_instance%coord_eff(:,o) + dr

                accepted = .TRUE.
            END IF
        END IF
    END SUBROUTINE mc_move

    SUBROUTINE mc_exchange(simulation_instance, energy_func)
        !========================================================================!
        ! This subroutine attempts to perform a grand canonical MC move by       !
        ! exchanging a particle with a reservoir                                 !
        !                                                                        !
        ! Reference:                                                             !
        ! Understanding Molecular Simulation                                     !
        ! Smit B and Frenkel D                                                   !
        ! Elsevier, 2002. https://doi.org/10.1016/B978-0-12-267351-1.X5000-7.    !
        ! pp. 32-33                                                              !
        !------------------------------------------------------------------------!
        ! simulation_instance      (inout) : simulation instance                 !
        ! energy_func                 (in) : pointer to energy function          !
        !========================================================================!
        IMPLICIT NONE
        TYPE(Simulation), INTENT(INOUT)      :: simulation_instance
        REAL(8)                              :: energy, prob_acceptance, arg
        REAL(8)                              :: r2, x2
        REAL(8)                              :: coord_tmp(3)
        INTEGER(4)                           :: o
        PROCEDURE(energy_function), POINTER  :: energy_func


        IF (RAND() .lt. 0.5) THEN
            IF (simulation_instance%npart .eq. 0) THEN
                ! There are no particles
                return
            END IF

            ! TODO: Test for overlap of particle upon insertion a displacement
            o = INT(simulation_instance%npart*RAND())+1

            ! Energy of particle to be removed
            energy = energy_func(simulation_instance, o, simulation_instance%coord(:,o))

            ! OLD ACCEPTANCE CRITERIA (Frenkel and Smith)
            !arg =  exp(energy/simulation_instance%temperature)
            !&* (simulation_instance%npart * 1d0) / (simulation_instance%pressure * simulation_instance%beta * simulation_instance%volume * 1d-30)

            ! New acceptance criteria (Tildesley)
            arg = energy/simulation_instance%temperature
            arg = arg + LOG( REAL(simulation_instance%npart) / (simulation_instance%activity*simulation_instance%volume*1d-30))

            IF (arg .gt. 0) THEN
                simulation_instance%coord(:,o) = simulation_instance%coord(:,simulation_instance%npart)
                simulation_instance%npart = simulation_instance%npart - 1
            ELSE IF (arg .gt. -30.0) THEN
                prob_acceptance = min(1.0d0,EXP(arg))
                IF ( RAND() .lt. prob_acceptance ) THEN
                    simulation_instance%coord(:,o) = simulation_instance%coord(:,simulation_instance%npart)
                    simulation_instance%npart = simulation_instance%npart - 1
                END IF
            END IF
        ELSE
            ! New particle at random position at the CNT
            r2 = simulation_instance%radius * simulation_instance%radius
            coord_tmp(1) = RAND()*2*(simulation_instance%radius-simulation_instance%sigma_cnt)-&
                    (simulation_instance%radius-simulation_instance%sigma_cnt)
            x2 = coord_tmp(1) * coord_tmp(1)
            coord_tmp(2) = RAND()*2*SQRT(r2-x2) - SQRT(r2 - x2)
            coord_tmp(3) = RAND()*simulation_instance%length

            ! Energy of particle to be inserted
            energy = energy_func(simulation_instance, simulation_instance%npart+1 , coord_tmp)

            ! New acceptance criteria (Tildesley)
            arg = - energy/simulation_instance%temperature
            arg = arg + LOG(simulation_instance%activity*simulation_instance%volume*1d-30/REAL(simulation_instance%npart+1))

            IF (arg .gt. 0) THEN
                simulation_instance%coord(:,simulation_instance%npart+1) = coord_tmp
                simulation_instance%npart = simulation_instance%npart + 1
            ELSE IF (arg .gt. -30.0) THEN
                prob_acceptance = min(1.0d0,EXP(arg))
                IF ( RAND() .lt. prob_acceptance ) THEN
                    simulation_instance%coord(:,simulation_instance%npart+1) = coord_tmp
                    simulation_instance%npart = simulation_instance%npart + 1
                END IF
            END IF
        END IF

    END SUBROUTINE mc_exchange
END MODULE monte_carlo_module
