!auto_modi:spll_albedo_from_nir_vis.D
MODULE MODI_ALBEDO_FROM_NIR_VIS
INTERFACE
      SUBROUTINE ALBEDO_FROM_NIR_VIS(PSW_BANDS,PALBNIR,PALBVIS,PALBUV,PDIR_ALB,PSCA_ALB)
IMPLICIT NONE
REAL, DIMENSION(:),   INTENT(IN) :: PSW_BANDS ! mean wavelength of each shortwave band (m)
REAL, DIMENSION(:),   INTENT(IN) :: PALBNIR   ! near infra-red albedo
REAL, DIMENSION(:),   INTENT(IN) :: PALBVIS   ! visible albedo
REAL, DIMENSION(:),   INTENT(IN) :: PALBUV    ! UV albedo
REAL, DIMENSION(:,:), INTENT(OUT):: PDIR_ALB  ! direct  albedo for each wavelength
REAL, DIMENSION(:,:), INTENT(OUT):: PSCA_ALB  ! diffuse albedo for each wavelength
END SUBROUTINE ALBEDO_FROM_NIR_VIS
END INTERFACE
END MODULE MODI_ALBEDO_FROM_NIR_VIS
