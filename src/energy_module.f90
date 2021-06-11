MODULE energy_module
    USE constants_module
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
        REAL(8)                :: lennard_jones_energy
        REAL(8)                :: rep_term, attract_term,sigma2

        ! Squared values
        sigma2 = sigma*sigma

        ! 4*epsilon((sigma/r)**12-(sigma/r)**6)
        attract_term = (sigma2 / r2) * (sigma2 / r2) * (sigma2 / r2)
        rep_term = attract_term * attract_term

        ! Lennard-Jones 12-6 energy
        lennard_jones_energy = 4.0d0 * epsilon * (rep_term - attract_term)

        RETURN
    END FUNCTION lennard_jones_energy

    FUNCTION lennard_jones_energy_tail_bulk_fluid(rcutoff, epsilon, sigma, density)
        !========================================================================!
        ! This function calculates the Lennard-Jones 12-6 energy tail correction !
        !                                                                        !
        ! Reference:                                                             !
        ! Understanding Molecular Simulation                                     !
        ! Smit B and Frenkel D                                                   !
        ! Elsevier, 2002. https://doi.org/10.1016/B978-0-12-267351-1.X5000-7.    !
        ! pp. 36-37                                                              !
        !------------------------------------------------------------------------!
        !------------------------------------------------------------------------!
        ! rcutoff      (in) : cutoff distance                                    !
        ! epsilon      (in) : potential well depth                               !
        ! sigma        (in) : distance at which the potential is zero            !
        ! density      (in) : density                                            !
        !========================================================================!
        REAL(8),    INTENT(IN) :: rcutoff, epsilon, sigma
        REAL(8)                :: lennard_jones_energy
        REAL(8)                :: rep_term, attract_term, sigma3 ,rcutoff3

        ! Squared values
        rcutoff3 = rcutoff*rcutoff
        sigma3 = sigma*sigma

        ! ((sigma/rcut)**9-(sigma/rcut)**3)
        attract_term = (sigma3 / rcutoff3)
        rep_term = (sigma3 / rcutoff3)*(sigma3 / rcutoff3)*(sigma3 / rcutoff3)

        ! Lennard-Jones 12-6 energy tail correction
        lennard_jones_energy = (8.0d0/3.0d0) * density * pi * epsilon * sigma3 * ((1.0d0/3.0d0)*rep_term - attract_term)

        RETURN
    END FUNCTION lennard_jones_energy_tail_bulk_fluid

    FUNCTION lennard_jones_energy_tail_fluid_solid(rcutoff, epsilon, sigma, density)
        !========================================================================!
        ! This function calculates the Lennard-Jones 12-6 energy tail correction !
        !                                                                        !
        ! Reference:                                                             !
        ! Siperstein, F.; Myers, A. L.; Talu, O.                                 !
        ! Long Range Corrections for Computer Simulations of Adsorption.         !
        ! Molecular Physics 2002, 100 (13), 2025–2030.                           !
        ! https://doi.org/10.1080/00268970110109916.                             !                                 !
        !------------------------------------------------------------------------!
        ! rcutoff      (in) : cutoff distance                                    !
        ! epsilon      (in) : potential well depth                               !
        ! sigma        (in) : distance at which the potential is zero            !
        ! density      (in) : density of solid atoms                             !
        !========================================================================!
        REAL(8),    INTENT(IN) :: rcutoff, epsilon, sigma
        REAL(8)                :: lennard_jones_energy
        REAL(8)                :: rep_term, attract_term, sigma3 ,rcutoff3

        ! Squared values
        rcutoff3 = rcutoff*rcutoff
        sigma3 = sigma*sigma

        ! ((sigma/rcut)**9-(sigma/rcut)**3)
        attract_term = (sigma3 / rcutoff3)
        rep_term = (sigma3 / rcutoff3)*(sigma3 / rcutoff3)*(sigma3 / rcutoff3)

        ! Lennard-Jones 12-6 energy tail correction
        ! Similar to bulk fluid correction but without the 1/2 factor that corrects for duplicate counting of
        ! pairwise interactioins in bulk fluids
        lennard_jones_energy = (16.0d0/3.0d0) * density * pi * epsilon * sigma3 * ((1.0d0/3.0d0)*rep_term - attract_term)

        RETURN
    END FUNCTION lennard_jones_energy_tail_fluid_solid

    FUNCTION total_system_energy(simulation_instance, energy_tail_correction)
        !========================================================================!
        ! This function calculates the total Lennard-Jones 12-6 energy           !
        ! for a atomistic nanotube representation                                !
        !                                                                        !
        ! energy_tail_correction is True if not present                          !
        !------------------------------------------------------------------------!
        ! simulation_instance    (in) : Simulation instance                      !
        ! energy_tail_correction (in) : flag to LJ energy tail correction        !
        !========================================================================!
        IMPLICIT NONE
        TYPE(Simulation),     INTENT(IN) :: simulation_instance
        LOGICAL(1), OPTIONAL, INTENT(IN) :: energy_tail_correction
        REAL(8)                          :: total_system_energy
        INTEGER(4)                       :: k

        IF (.not. PRESENT(energy_tail_correction)) THEN
            energy_tail_correction = .true.
        END IF

        ! Total system energy is equal to sum of invidual energy particle contributions divided by 2
        total_system_energy = 0
        DO k=1,simulation_instance%npart
            total_system_energy = total_system_energy &
                    + particle_energy_atomistic(simulation_instance, k, simulation_instance%coord(:,k)) / 2.0d0
        END DO


        ! Add Lennard-Jones 12-6 Tail correction
        IF (energy_tail_correction .eqv. .True.) THEN
            ! Add bulk fluid correction
            total_system_energy = total_system_energy &
                    + simulation_instance%npart &
                    * lennard_jones_energy_tail_fluid_solid(simulation_instance%fluid_cut_off, &
                                    simulation_instance%eps_fluid, simulation_instance%sigma_fluid,simulation_instance%density_fluid)
            ! Add fluid-solid correction
            total_system_energy = total_system_energy &
                    + simulation_instance%ncnt &
                            * lennard_jones_energy_tail_fluid_solid(simulation_instance%fluid_cut_off, &
                                    simulation_instance%eps_cnt, simulation_instance%sigma_cnt,simulation_instance%density_cnt)
        END IF

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
            CALL pbc_z_distance(r, simulation_instance%length)
            r2 = SUM(r*r)

            ! If within cut-off take interaction into account
            IF (r2 .lt. simulation_instance%fluid_cut_off * simulation_instance%fluid_cut_off) THEN
                energy_particle_fluid = energy_particle_fluid + &
                        lennard_jones_energy(r2, simulation_instance%eps_fluid, simulation_instance%sigma_fluid )
            END IF
        END DO

        DO k=index+1,simulation_instance%npart
            r = coord(:)-simulation_instance%coord(:,k)
            CALL pbc_z_distance(r, simulation_instance%length)

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
            CALL pbc_z_distance(r, simulation_instance%length)

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
            CALL pbc_z_distance(r, simulation_instance%length)

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
            CALL pbc_z_distance(r, simulation_instance%length)

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
