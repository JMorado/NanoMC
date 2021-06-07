MODULE energy_module
  USE pbc_displacement_module
  USE simulation_module
  IMPLICIT NONE



CONTAINS
  FUNCTION lj_fluid(this,r2)
    !----------------------------------------------------------------------------
    !TODO:make it work with/without reduced units
    !
    !----------------------------------------------------------------------------
    CLASS(Simulation), INTENT(IN) :: this

    REAL(8), INTENT(IN) :: r2
    REAL(8)             :: lj_fluid, attract_term, rep_term, sigma2

    sigma2 = this%sigma_fluid * this%sigma_fluid                               ! Si
    attract_term = (sigma2 / r2) * (sigma2 / r2) * (sigma2 / r2)
    rep_term = attract_term * attract_term

    lj_fluid = 4.0d0 * this%eps_fluid * (rep_term - attract_term)

    RETURN
  END FUNCTION lj_fluid

  FUNCTION lj_fluid_cnt(this, r2)
    !----------------------------------------------------------------------------
    !
    !
    !----------------------------------------------------------------------------
    CLASS(Simulation), INTENT(IN) :: this
    REAL(8), INTENT(IN) :: r2
    REAL(8)             :: lj_fluid_cnt, attract_term, rep_term, sigma2

    sigma2 = this%sigma_cnt * this%sigma_cnt                               ! Sigma squared
    attract_term = (sigma2 / r2) * (sigma2 / r2) * (sigma2 / r2)
    rep_term = attract_term * attract_term

    lj_fluid_cnt = 4.0d0 * this%eps_cnt * (rep_term - attract_term)

    RETURN
  END FUNCTION lj_fluid_cnt


  FUNCTION lj_fluid_cnt_atomistic(this, r2)
    !----------------------------------------------------------------------------
    !
    !
    !----------------------------------------------------------------------------
    CLASS(Simulation), INTENT(IN) :: this
    REAL(8), INTENT(IN) :: r2
    REAL(8)             :: lj_fluid_cnt_atomistic, attract_term, rep_term, sigma2

    sigma2 = this%sigma_cnt * this%sigma_cnt                               ! Sigma squared
    attract_term = (sigma2 / r2) * (sigma2 / r2) * (sigma2 / r2)
    rep_term = attract_term * attract_term

    lj_fluid_cnt_atomistic = 4.0d0 * this%eps_cnt * (rep_term - attract_term)

    RETURN
  END FUNCTION lj_fluid_cnt_atomistic


  FUNCTION system_energy(this, coord, sigma, eps)
    IMPLICIT NONE
    !----------------------------------------------------------------------------
    ! TODO: check this subroutine
    !
    !----------------------------------------------------------------------------
    CLASS(Simulation), INTENT(IN) :: this
    REAL(8),           INTENT(IN) :: sigma, eps
    REAL(8),           INTENT(IN) :: coord(:)
    REAL(8)                       :: system_energy, energy_cnt_fluid, energy_fluid, r2
    INTEGER(4)                    :: k, l


    system_energy = 0
    DO k=1,this%npart
       energy_fluid = particle_energy(this, k, coord) / 2.0d0

       r2 = this%radius*this%radius - SUM( coord(1:2)*coord(1:2) )
       energy_cnt_fluid = lj_fluid_cnt(this,r2)

       system_energy = system_energy + energy_fluid + energy_cnt_fluid
    END DO

    RETURN
  END FUNCTION system_energy


  FUNCTION particle_energy(this, index, coord)
    IMPLICIT NONE
    !-----------------------------------------------------------------------
    ! This functions returns the energy of a particle with index 'index'
    ! at the given postion 'coord' in the CNT continuum model.
    !-----------------------------------------------------------------------
    CLASS(Simulation), INTENT(IN)      :: this
    REAL(8),           INTENT(IN)      :: coord(:)
    INTEGER(4),        INTENT(IN)      :: index
    REAL(8)                            :: r(3)
    REAL(8)                            :: particle_energy, energy_cnt_part, energy_part, r2
    INTEGER(4)                         :: k


    particle_energy = 0
    ! Fluid particle-particle contribution to the energy
    DO k=1,index-1
       r = coord(:)-this%coord(:,k)
       CALL pbc_z(r, this%length)

       r2 = SUM(r*r)
       particle_energy = particle_energy + lj_fluid(this,r2)
    END DO

    DO k=index+1,this%npart
       r = coord(:)-this%coord(:,k)

       CALL pbc_z(r, this%length)
       r2 = SUM(r*r)

       particle_energy = particle_energy + lj_fluid(this,r2)
    END DO


    ! Fluid-CNT contribution to the energy
    r2 = this%radius - SQRT(SUM( coord(1:2)*coord(1:2) ))
    r2 = r2*r2
    energy_cnt_part = lj_fluid_cnt(this,r2)

    ! Add up different contributions to the energy
    ! E_particle_total = E_particle_particle + E_particle_cnt
    particle_energy = particle_energy + energy_cnt_part

    RETURN
  END FUNCTION PARTICLE_ENERGY



  FUNCTION particle_energy_atomistic(this, index, coord)
    IMPLICIT NONE
    !-----------------------------------------------------------------------
    ! This functions returns the energy of a particle with index 'index'
    ! at the given postion 'coord' in the CNT atomistic model.
    !-----------------------------------------------------------------------
    
    CLASS(Simulation), INTENT(IN) :: this
    REAL(8),           INTENT(IN) :: coord(:)
    INTEGER(4),        INTENT(IN) :: index
    REAL(8)                       :: r(3)
    REAL(8)                       :: particle_energy_atomistic, energy_cnt_part, energy_part, r2
    INTEGER(4)                    :: k,l
   

    particle_energy_atomistic = 0
    ! Fluid inter-particle contribution to the energy
    DO k=1,index-1
       ! Compute distance r
       r = coord(:)-this%coord(:,k)
       ! Apply PBC to r
       CALL pbc_z(r, this%length)
       ! Compute r squared of the periodic distance
       r2 = SUM(r*r)

       particle_energy_atomistic = particle_energy_atomistic + lj_fluid(this,r2)
    END DO
 
    DO k=index+1,this%npart
       ! Compute distance r
       r = coord(:)-this%coord(:,k)
       ! Apply PBC to r
       CALL pbc_z(r, this%length)
       ! Compute r squared of the periodic distance
       r2 = SUM(r*r)

       particle_energy_atomistic = particle_energy_atomistic + lj_fluid(this,r2)
    END DO

  
    ! Fluid-CNT contribution to the energy
    DO k=1,this%ncnt
       ! Compute distance r       
       r = coord(:)-this%coord_cnt(:,k)
       ! Apply PBC to r
       CALL pbc_z(r, this%length)
       ! Compute r squared of the periodic distance
       r2 = SUM(r*r)

       particle_energy_atomistic = particle_energy_atomistic + lj_fluid_cnt_atomistic(this,r2)
    END DO

    RETURN
  END FUNCTION PARTICLE_ENERGY_ATOMISTIC

  FUNCTION virial(this)
    !TODO: does not seem to be working properly yet unfortunately
    IMPLICIT NONE
    ! Input class
    CLASS(Simulation), INTENT(IN) :: this
    ! Intermediate variables
    INTEGER(4) :: k, l
    REAL(8) :: r(3), r2
    REAL(8) :: attract_term, rep_term
    REAL(8) :: sigmaf2, sigmacnt2
    ! Output variable
    REAL(8) :: virial


    sigmaf2 = this%sigma_fluid*this%sigma_fluid          ! Sigma fluid squared
    sigmacnt2 = this%sigma_cnt*this%sigma_cnt            ! Sigma CNT squared


    virial = 0
    DO k=1,this%npart
       DO l=k+1,this%npart
          r = this%coord(:,k)-this%coord(:,l)
          r(3) = r(3)-this%length*NINT(r(3)/this%length)
          r2 = SUM(r*r)

          attract_term = (sigmaf2 / r2) * (sigmaf2 / r2) * (sigmaf2 / r2)
          rep_term = attract_term * attract_term


          virial = virial + 48 * this%eps_fluid * (rep_term - 0.5 * attract_term) / SQRT(r2)
       END DO
       ! Fluid-CNT contribution to the eergy
       r2 = this%radius - SQRT(SUM( this%coord(1:2,k)*this%coord(1:2,k) ))
       r2 = r2*r2
       attract_term = (sigmacnt2 / r2) * (sigmacnt2 / r2) * (sigmacnt2 / r2)
       rep_term = attract_term * attract_term

      
       virial = virial + 48 * this%eps_cnt * (rep_term - 0.5 * attract_term) / SQRT(r2)
    END DO

    virial = (1.0d0/3.0d0) * virial / REAL(this%npart)


    RETURN
  END FUNCTION virial
END MODULE ENERGY_MODULE
