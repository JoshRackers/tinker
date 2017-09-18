c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module potderivs  --  electric potential, field, field  ##        
c     ##                        gradient and field hessian        ##
c     ##                                                          ##
c     ##############################################################
c
c
c     potm        electric potential at every multipole site
c     fieldm      electric field with mscale exclusion rules
c     fieldd      electric field with dscale exclusion rules
c     fieldp      electric field with pscale exclusion rules
c     gradfieldm  electric field gradient with mscale exclusion rules
c     gradfieldp  electric field gradient with pscale exclusion rules
c     hessfieldm  electric field hession with mscale exclusion rules
c
c
      module potderivs
      implicit none
      real*8, allocatable :: pot(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: gradfield(:,:,:)
      real*8, allocatable :: hessfield(:,:,:,:)
c
      real*8, allocatable :: potm(:)
      real*8, allocatable :: fieldm(:,:)
      real*8, allocatable :: gradfieldm(:,:,:)
      real*8, allocatable :: hessfieldm(:,:,:,:)
c
      real*8, allocatable :: pot_ewald(:)
      real*8, allocatable :: field_ewald(:,:)
      real*8, allocatable :: gradfield_ewald(:,:,:)
      real*8, allocatable :: hessfield_ewald(:,:,:,:)
c
      real*8, allocatable :: fieldd_thole(:,:)
      real*8, allocatable :: fieldp_thole(:,:)
      real*8, allocatable :: gradfieldd_thole(:,:,:)
      real*8, allocatable :: gradfieldp_thole(:,:,:)
c
      real*8, allocatable :: fieldd_func(:,:)
      real*8, allocatable :: fieldp_func(:,:)
      real*8, allocatable :: gradfieldd_func(:,:,:)
      real*8, allocatable :: gradfieldp_func(:,:,:)
c
      real*8, allocatable :: potm_gordon(:)
      real*8, allocatable :: fieldm_gordon(:,:)
      real*8, allocatable :: gradfieldm_gordon(:,:,:)
      real*8, allocatable :: hessfieldm_gordon(:,:,:,:)
c
      real*8, allocatable :: potm_goverlap(:)
      real*8, allocatable :: fieldm_goverlap(:,:)
      real*8, allocatable :: gradfieldm_goverlap(:,:,:)
      real*8, allocatable :: hessfieldm_goverlap(:,:,:,:)
c
      real*8, allocatable :: nucpotm_gordon(:)
      real*8, allocatable :: nucfieldm_gordon(:,:)
c
      real*8, allocatable :: fieldd_gordon(:,:)
      real*8, allocatable :: fieldp_gordon(:,:)
      real*8, allocatable :: gradfieldd_gordon(:,:,:)
      real*8, allocatable :: gradfieldp_gordon(:,:,:)
c
      real*8, allocatable :: udnucfieldp_gordon(:,:)
      real*8, allocatable :: upnucfieldd_gordon(:,:)
c
      real*8, allocatable :: udfield_gordon(:,:)
      real*8, allocatable :: upfield_gordon(:,:)
c
      real*8, allocatable :: udfieldmu_gordon(:,:)
      real*8, allocatable :: upfieldmu_gordon(:,:)
c
      real*8, allocatable :: udgradfield_gordon(:,:,:)
      real*8, allocatable :: upgradfield_gordon(:,:,:)
      real*8, allocatable :: udhessfield_gordon(:,:,:,:)
      real*8, allocatable :: uphessfield_gordon(:,:,:,:)
c
      real*8, allocatable :: udgradfieldmu_gordon(:,:,:)
      real*8, allocatable :: upgradfieldmu_gordon(:,:,:)
      real*8, allocatable :: udhessfieldmu_gordon(:,:,:,:)
      real*8, allocatable :: uphessfieldmu_gordon(:,:,:,:)
c
      real*8, allocatable :: udfieldp_gordon(:,:)
      real*8, allocatable :: udgradfieldp_gordon(:,:,:)
      real*8, allocatable :: udhessfieldp_gordon(:,:,:,:)
c
      real*8, allocatable :: upfieldd_gordon(:,:)
      real*8, allocatable :: upgradfieldd_gordon(:,:,:)
      real*8, allocatable :: uphessfieldd_gordon(:,:,:,:)
c
      real*8, allocatable :: fieldd_gordonreg(:,:)
      real*8, allocatable :: fieldp_gordonreg(:,:)
      real*8, allocatable :: gradfieldd_gordonreg(:,:,:)
c
      real*8, allocatable :: potm_piquemal(:)
      real*8, allocatable :: fieldm_piquemal(:,:)
      real*8, allocatable :: gradfieldm_piquemal(:,:,:)
      real*8, allocatable :: hessfieldm_piquemal(:,:,:,:)
c
      real*8, allocatable :: nucpotm_piquemal(:)
      real*8, allocatable :: nucfieldm_piquemal(:,:)
c
      real*8, allocatable :: fieldd_piquemal(:,:)
      real*8, allocatable :: fieldp_piquemal(:,:)
      real*8, allocatable :: gradfieldd_piquemal(:,:,:)
      real*8, allocatable :: gradfieldp_piquemal(:,:,:)
c
      real*8, allocatable :: udnucfieldp_piquemal(:,:)
      real*8, allocatable :: upnucfieldd_piquemal(:,:)
c
      real*8, allocatable :: udfield_piquemal(:,:)
      real*8, allocatable :: upfield_piquemal(:,:)
      real*8, allocatable :: udgradfield_piquemal(:,:,:)
      real*8, allocatable :: upgradfield_piquemal(:,:,:)
      real*8, allocatable :: udhessfield_piquemal(:,:,:,:)
      real*8, allocatable :: uphessfield_piquemal(:,:,:,:)
c
      real*8, allocatable :: udfieldp_piquemal(:,:)
      real*8, allocatable :: udgradfieldp_piquemal(:,:,:)
      real*8, allocatable :: udhessfieldp_piquemal(:,:,:,:)
c
      real*8, allocatable :: upfieldd_piquemal(:,:)
      real*8, allocatable :: upgradfieldd_piquemal(:,:,:)
      real*8, allocatable :: uphessfieldd_piquemal(:,:,:,:)
c
      real*8, allocatable :: fieldd_piquemalreg(:,:)
      real*8, allocatable :: fieldp_piquemalreg(:,:)
c
c      real*8, allocatable :: potm_piquemalone(:)
c      real*8, allocatable :: fieldm_piquemalone(:,:)
c      real*8, allocatable :: fieldd_piquemalone(:,:)
c      real*8, allocatable :: fieldp_piquemalone(:,:)
c      real*8, allocatable :: gradfieldm_piquemalone(:,:,:)
c      real*8, allocatable :: gradfieldp_piquemalone(:,:,:)
c      real*8, allocatable :: potm_piquemaltwo(:)
c      real*8, allocatable :: fieldm_piquemaltwo(:,:)
c      real*8, allocatable :: fieldd_piquemaltwo(:,:)
c      real*8, allocatable :: fieldp_piquemaltwo(:,:)
c      real*8, allocatable :: gradfieldm_piquemaltwo(:,:,:)
c      real*8, allocatable :: gradfieldp_piquemaltwo(:,:,:)
cccccccccccccccccccc
c
c     all induced dipole field quantities are d or p (coming from
c     d or p dipoles)
c
      real*8, allocatable :: udfield(:,:)
      real*8, allocatable :: upfield(:,:)
      real*8, allocatable :: udgradfield(:,:,:)
      real*8, allocatable :: upgradfield(:,:,:)
      real*8, allocatable :: udhessfield(:,:,:,:)
      real*8, allocatable :: uphessfield(:,:,:,:)
c
      real*8, allocatable :: udfield_ewald(:,:)
      real*8, allocatable :: upfield_ewald(:,:)
      real*8, allocatable :: udgradfield_ewald(:,:,:)
      real*8, allocatable :: upgradfield_ewald(:,:,:)
      real*8, allocatable :: udhessfield_ewald(:,:,:,:)
      real*8, allocatable :: uphessfield_ewald(:,:,:,:)
c
      real*8, allocatable :: udfield_thole(:,:)
      real*8, allocatable :: upfield_thole(:,:)
      real*8, allocatable :: udgradfield_thole(:,:,:)
      real*8, allocatable :: upgradfield_thole(:,:,:)
      real*8, allocatable :: udhessfield_thole(:,:,:,:)
      real*8, allocatable :: uphessfield_thole(:,:,:,:)
c
      real*8, allocatable :: upfieldd_thole(:,:)
      real*8, allocatable :: upgradfieldd_thole(:,:,:)
      real*8, allocatable :: uphessfieldd_thole(:,:,:,:)
c
      real*8, allocatable :: udfieldp_thole(:,:)
      real*8, allocatable :: udgradfieldp_thole(:,:,:)
      real*8, allocatable :: udhessfieldp_thole(:,:,:,:)
c
      real*8, allocatable :: udfield_func(:,:)
      real*8, allocatable :: upfield_func(:,:)
      real*8, allocatable :: udgradfield_func(:,:,:)
      real*8, allocatable :: upgradfield_func(:,:,:)
      real*8, allocatable :: udhessfield_func(:,:,:,:)
      real*8, allocatable :: uphessfield_func(:,:,:,:)
c
      real*8, allocatable :: upfieldd_func(:,:)
      real*8, allocatable :: upgradfieldd_func(:,:,:)
      real*8, allocatable :: uphessfieldd_func(:,:,:,:)
c     
      real*8, allocatable :: udfieldp_func(:,:)
      real*8, allocatable :: udgradfieldp_func(:,:,:)
      real*8, allocatable :: udhessfieldp_func(:,:,:,:)
c
c     reciprocal space field quantities
c
      real*8, allocatable :: pot_recip(:)
      real*8, allocatable :: field_recip(:,:)
      real*8, allocatable :: gradfield_recip(:,:,:)
      real*8, allocatable :: hessfield_recip(:,:,:,:)
c
c     the d and p here only refer to the two types of induced dipoles
c     not applying the d and p scale factors
c
      real*8, allocatable :: udfield_recip(:,:)
      real*8, allocatable :: upfield_recip(:,:)
      real*8, allocatable :: udgradfield_recip(:,:,:)
      real*8, allocatable :: upgradfield_recip(:,:,:)
      real*8, allocatable :: udhessfield_recip(:,:,:,:)
      real*8, allocatable :: uphessfield_recip(:,:,:,:)
      real*8, allocatable :: udphessfield_recip(:,:,:,:)
c
c     self fields
c
      real*8, allocatable :: field_self(:,:)
      real*8, allocatable :: udfield_self(:,:)
      real*8, allocatable :: upfield_self(:,:)
c
c     logical flags for types of damping needed
c
      logical damp_none
      logical damp_ewald
      logical damp_thole
      logical damp_func
      logical damp_gordon
      logical damp_piquemal
      logical damp_gordonreg
      logical damp_goverlap
      save
      end
