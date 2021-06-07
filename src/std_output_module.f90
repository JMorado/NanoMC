MODULE std_output_module
  USE simulation_module
  IMPLICIT NONE



CONTAINS
  SUBROUTINE STD_OUTPUT_INITIALIZE()
    IMPLICIT NONE

    WRITE(6,*) "!----------------------------------------------------------------------------!"
    WRITE(6,*) "!                   _   _                   __  __  ____                     !"
    WRITE(6,*) "!                  | \ | |                 |  \/  |/ ____|                   !"
    WRITE(6,*) "!                  |  \| | __ _ _ __   ___ | \  / | |                        !"
    WRITE(6,*) "!                  | . ` |/ _` | '_ \ / _ \| |\/| | |                        !"
    WRITE(6,*) "!                  | |\  | (_| | (_| | (_) | |  | | |____                    !"
    WRITE(6,*) "!                  |_| \_|\__,_|_| |_|\___/|_|  |_|\_____|                   !"
    WRITE(6,*) "!                                                                            !"
    WRITE(6,*) "!----------------------------------------------------------------------------!"
    WRITE(6,*) "!                        NanoMC code by Joao Morado                          !"
    WRITE(6,*) "!----------------------------------------------------------------------------!"
  END SUBROUTINE STD_OUTPUT_INITIALIZE


  SUBROUTINE STD_OUTPUT_SIMULATION_DETAILS(this)
    IMPLICIT NONE
    CLASS(Simulation), INTENT(IN) :: this

    WRITE(6,*) "!----------------------------------------------------------------------------!"
    WRITE(6,*) "!                            Simulation Details                              !"
    WRITE(6,*) "!----------------------------------------------------------------------------!"

    WRITE(6,*)
    WRITE(6,'(A)') " * Simulation Variables *"
    WRITE(6,'(30A,A)') "   - CNT model: ", this%model
    WRITE(6,'(30A,A)') "   - Ensemble: ",this%ensemble

    WRITE(6,*)
    WRITE(6,'(A)') " * CNT Parameters *"
    WRITE(6,'(A30,F8.4)') adjustl("   - CNT length (A): "),this%length
    WRITE(6,'(A30,F8.4)') adjustl("   - CNT radius (A): "),this%radius
    WRITE(6,'(A)') " * "

    WRITE(6,*)
    WRITE(6,'(A)') " * Lennard-Jones 12-6 parameters *"
    WRITE(6,'(A30,F8.4)') adjustl("   - fluid-fluid cut-off (A): "),this%fluid_cut_off
    WRITE(6,'(A30,F8.4)') adjustl("   - CNT-CNTT cut-off (A): "),this%cnt_cut_off
    WRITE(6,'(A30,F8.4)') adjustl("   - Fluid sigma (A): "),this%sigma_fluid
    WRITE(6,'(A30,F8.4)') adjustl("   - Fluid epsilon: "),this%eps_fluid
    WRITE(6,'(A30,F8.4)') adjustl("   - CNT sigma (A): "),this%sigma_cnt
    WRITE(6,'(A30,F8.4)') adjustl("   - CNT epsilon: "),this%eps_cnt
    WRITE(6,'(A)') " * "

    WRITE(6,*)
    WRITE(6,'(A)') " * Seed variables *"
    WRITE(6,'(A30,A)') adjustl("   - Output seed: "),this%output_seed
    WRITE(6,'(A30,I8.4)') adjustl("   - Random seed: "),this%seed
    WRITE(6,'(A)') " * "

    WRITE(6,'(A)') " * Output variables *"
    WRITE(6,'(A30,I8.4)') adjustl("   - ntwx: "), this%ntwx
    WRITE(6,'(A30,I8.4)') adjustl("   - ntpr: "), this%ntpr
    WRITE(6,'(A)') " * "

    WRITE(6,'(A)') " * Input configuration variables"
    WRITE(6,'(A30,A)') adjustl("   - Initial configuration: "), this%initial_config_type
    WRITE(6,'(A)') " * "

    WRITE(6,*) "!----------------------------------------------------------------------------!"
    WRITE(6,*) "!                         End of Simulation Details                          !"
    WRITE(6,*) "!----------------------------------------------------------------------------!"
    WRITE(6,*)

  END SUBROUTINE STD_OUTPUT_SIMULATION_DETAILS

  SUBROUTINE STD_OUTPUT_STARTING_SIM()
    IMPLICIT NONE


    WRITE(6,*) "!----------------------------------------------------------------------------!"
    WRITE(6,*) "!                                                                            !"
    WRITE(6,*) "!                  Starting the uVT ensemble simulation                      !"
    WRITE(6,*) "!                                                                            !"
    WRITE(6,*) "!----------------------------------------------------------------------------!"

  END SUBROUTINE STD_OUTPUT_STARTING_SIM


  
  SUBROUTINE STD_OUTPUT_END_SIM()
    IMPLICIT NONE

    WRITE(6,*) "!----------------------------------------------------------------------------!"
    WRITE(6,*) "!                                                                            !"
    WRITE(6,*) "!                     End of uVT ensemble simulation                         !"
    WRITE(6,*) "!                                                                            !"
    WRITE(6,*) "!----------------------------------------------------------------------------!"


  END SUBROUTINE STD_OUTPUT_END_SIM





END MODULE std_output_module
