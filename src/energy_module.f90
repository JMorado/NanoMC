MODULE energy_module
    USE pbc_module
    USE simulation_module
    IMPLICIT NONE

CONTAINS
    FUNCTION lennard_jones_energy(r2, epsilon, sigma)
        !========================================================================!
        ! This function calculates the Lennard-Jones 12-6 energy                 !
        !------------------------------------------------------------------------!
        ! r2           (in) : distance squared                                   !
        ! epsilon      (in) : potential well depth                               !
        ! sigma        (in) : distance at which the potential is zero            !
        !========================================================================!
        REAL(8),    INTENT(IN) :: r2, epsilon, sigma
        REAL(8)                :: lennard_jones_energy, rep_term, attract_term,sigma2

        ! Squared values
        sigma2 = sigma*sigma

        ! 4*epsilon((sigma/r)**12-(sigma/r)**6)
        attract_term = (sigma2 / r2) * (sigma2 / r2) * (sigma2 / r2)
        rep_term = attract_term * attract_term

        ! Lennard-Jones 12-6 energy
        lennard_jones_energy = 4.0d0 * epsilon * (rep_term - attract_term)

        RETURN
    END FUNCTION lennard_jones_energy

    FUNCTION total_system_energy(simulation_instance)
        !========================================================================!
        ! This function calculates the total Lennard-Jones 12-6 energy           !
        ! for a atomistic nanotube representation                                !
        !------------------------------------------------------------------------!
        ! simulation_instance (in) : Simulation instance                         !
        !========================================================================!
        IMPLICIT NONE
        TYPE(Simulation),  INTENT(IN) :: simulation_instance
        REAL(8)                       :: total_system_energy
        INTEGER(4)                    :: k

        ! Total system energy is equal to sum of invidual energy particle contributions divided by 2
        total_system_energy = 0
        DO k=1,simulation_instance%npart
            total_system_energy = total_system_energy &
                    + particle_energy_atomistic(simulation_instance, k, simulation_instance%coord(:,k)) / 2.0d0
        END DO

        RETURN
    END FUNCTION total_system_energy

    FUNCTION particle_energy_continuum(simulation_instance, index, coord)
        IMPLICIT NONE
        !========================================================================!
        ! This function calculates the particle Lennard-Jones 12-6 energy        !
        ! for a continuum nanotube representation                                !
        !------------------------------------------------------------------------!
        ! simulation_instance (in) : Simulation instance                         !
        ! index               (in) : index of the particle                       !
        ! coord(3)            (in) : coordinates of the particle                 !
        !========================================================================!
        TYPE(Simulation),  INTENT(IN)      :: simulation_instance
        REAL(8),           INTENT(IN)      :: coord(:)
        INTEGER(4),        INTENT(IN)      :: index
        REAL(8)                            :: r(3)
        REAL(8)                            :: particle_energy_continuum, energy_nt_particle, energy_particle_fluid, r2
        INTEGER(4)                         :: k

        !========================================================================!
        !                         Fluid-fluid LJ 12-6 energy                     !
        !========================================================================!
        energy_particle_fluid = 0

        DO k=1,index-1
            r = coord(:)-simulation_instance%coord(:,k)
            CALL pbc(r, (/simulation_instance%cell_instance%a,simulation_instance%cell_instance%b,&
                    simulation_instance%cell_instance%c/))
            r2 = SUM(r*r)

            ! If within cut-off take interaction into account
            IF (r2 .lt. simulation_instance%fluid_cut_off * simulation_instance%fluid_cut_off) THEN
                energy_particle_fluid = energy_particle_fluid + &
                        lennard_jones_energy(r2, simulation_instance%eps_fluid, simulation_instance%sigma_fluid )
            END IF
        END DO

        DO k=index+1,simulation_instance%npart
            r = coord(:)-simulation_instance%coord(:,k)
            CALL pbc(r, (/simulation_instance%cell_instance%a,simulation_instance%cell_instance%b,&
                    simulation_instance%cell_instance%c/))
            r2 = SUM(r*r)

            ! If within cut-off take interaction into account
            IF (r2 .lt. simulation_instance%fluid_cut_off * simulation_instance%fluid_cut_off) THEN
                energy_particle_fluid = energy_particle_fluid + &
                        lennard_jones_energy(r2, simulation_instance%eps_fluid, simulation_instance%sigma_fluid )
            END IF
        END DO

        !========================================================================!
        !                      Fluid-NT(continuum) LJ 12-6 energy                !
        !========================================================================!
        r2 = simulation_instance%radius - SQRT(SUM( coord(1:2)*coord(1:2) ))
        r2 = r2*r2

        IF (r2 .lt. simulation_instance%cnt_cut_off * simulation_instance%cnt_cut_off) THEN
            energy_nt_particle = energy_nt_particle &
                    + lennard_jones_energy(r2, simulation_instance%eps_cnt, simulation_instance%sigma_cnt)
        END IF

        !========================================================================!
        !                            Total LJ 12-6 energy                        !
        !========================================================================!
        ! Add up different contributions to the energy
        ! E_particle_total = E_particle_particle + E_particle_cnt
        particle_energy_continuum = energy_particle_fluid + energy_nt_particle

        RETURN
    END FUNCTION particle_energy_continuum

    FUNCTION particle_energy_atomistic(simulation_instance, index, coord)
        !========================================================================!
        ! This function calculates the particle Lennard-Jones 12-6 energy        !
        ! for an atomistic nanotube representation                               !
        !------------------------------------------------------------------------!
        ! simulation_instance (in) : Simulation instance                         !
        ! index               (in) : index of the particle                       !
        ! coord(3)            (in) : coordinates of the particle                 !
        !========================================================================!
        IMPLICIT NONE
        TYPE(Simulation),  INTENT(IN) :: simulation_instance
        REAL(8),           INTENT(IN) :: coord(:)
        INTEGER(4),        INTENT(IN) :: index
        REAL(8)                       :: r(3)
        REAL(8)                       :: particle_energy_atomistic, energy_nt_particle, energy_particle_fluid, r2
        INTEGER(4)                    :: k, l

        !========================================================================!
        !                         Fluid-fluid LJ 12-6 energy                     !
        !========================================================================!
        energy_particle_fluid = 0

        ! Fluid inter-particle contribution to the energy
        DO k=1,index-1
            ! Compute distance r
            r = coord(:)-simulation_instance%coord(:,k)

            ! Apply PBC
            CALL pbc(r, (/simulation_instance%cell_instance%a,simulation_instance%cell_instance%b,&
                    simulation_instance%cell_instance%c/))

            ! Compute r squared of the periodic distance
            r2 = SUM(r*r)

            ! If within cut-off radius take interaction into account
            IF (r2 .lt. simulation_instance%fluid_cut_off * simulation_instance%fluid_cut_off) THEN
                energy_particle_fluid = energy_particle_fluid + &
                        lennard_jones_energy(r2, simulation_instance%eps_fluid, simulation_instance%sigma_fluid )
            END IF
        END DO

        DO k=index+1,simulation_instance%npart
            ! Compute distance r
            r = coord(:)-simulation_instance%coord(:,k)
            ! Apply PBC to r
            CALL pbc(r, (/simulation_instance%cell_instance%a,simulation_instance%cell_instance%b,&
                    simulation_instance%cell_instance%c/))
            ! Compute r squared of the periodic distance
            r2 = SUM(r*r)

            ! If within cut-off radius take interaction into account
            IF (r2 .lt. simulation_instance%fluid_cut_off * simulation_instance%fluid_cut_off) THEN
                energy_particle_fluid = energy_particle_fluid + &
                        lennard_jones_energy(r2, simulation_instance%eps_fluid, simulation_instance%sigma_fluid )
            END IF
        END DO

        !========================================================================!
        !                      Fluid-NT(atomistic) LJ 12-6 energy                !
        !========================================================================!
        energy_nt_particle = 0
        DO k=1,simulation_instance%ncnt
            ! Compute distance r
            r = coord(:)-simulation_instance%coord_cnt(:,k)

            ! Apply PBC
            CALL pbc(r, (/simulation_instance%cell_instance%a,simulation_instance%cell_instance%b,&
                    simulation_instance%cell_instance%c/))

            ! Compute r squared of the periodic distance
            r2 = SUM(r*r)

            IF (r2 .lt. simulation_instance%cnt_cut_off * simulation_instance%cnt_cut_off) THEN
                energy_nt_particle = energy_nt_particle &
                        + lennard_jones_energy(r2, simulation_instance%eps_cnt, simulation_instance%sigma_cnt)
            END IF
        END DO

        !========================================================================!
        !                            Total LJ 12-6 energy                        !
        !========================================================================!
        particle_energy_atomistic = energy_nt_particle + energy_particle_fluid

        RETURN
    END FUNCTION particle_energy_atomistic
END MODULE energy_module
