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
c     clean these up later
c
c      real*8, allocatable :: potm_gordonone(:)
c      real*8, allocatable :: fieldm_gordonone(:,:)
c      real*8, allocatable :: fieldd_gordonone(:,:)
c      real*8, allocatable :: fieldp_gordonone(:,:)
c      real*8, allocatable :: gradfieldm_gordonone(:,:,:)
c      real*8, allocatable :: gradfieldp_gordonone(:,:,:)
c      real*8, allocatable :: potm_gordontwo(:)
c      real*8, allocatable :: fieldm_gordontwo(:,:)
c      real*8, allocatable :: fieldd_gordontwo(:,:)
c      real*8, allocatable :: fieldp_gordontwo(:,:)
c      real*8, allocatable :: gradfieldm_gordontwo(:,:,:)
c      real*8, allocatable :: gradfieldp_gordontwo(:,:,:)
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
      real*8, allocatable :: upfieldp_thole(:,:)
      real*8, allocatable :: upgradfieldp_thole(:,:,:)
      real*8, allocatable :: uphessfieldp_thole(:,:,:,:)
c
      real*8, allocatable :: udfieldd_thole(:,:)
      real*8, allocatable :: udgradfieldd_thole(:,:,:)
      real*8, allocatable :: udhessfieldd_thole(:,:,:,:)
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
      logical damp_gordonone
      logical damp_gordontwo
      logical damp_piquemalone
      logical damp_piquemaltwo
      save
      end
