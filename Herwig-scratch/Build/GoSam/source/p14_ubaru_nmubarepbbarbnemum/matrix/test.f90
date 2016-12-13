program test
   use p14_ubaru_nmubarepbbarbnemum_config, only: ki, logfile, nlo_prefactors
   use p14_ubaru_nmubarepbbarbnemum_kinematics, only: dotproduct, boost_to_cms
   use p14_ubaru_nmubarepbbarbnemum_model, only: parse
   use p14_ubaru_nmubarepbbarbnemum_matrix, only: samplitude, &
     & initgolem, exitgolem, ir_subtraction
   use p14_ubaru_nmubarepbbarbnemum_color, only: numcs, CA
   use p14_ubaru_nmubarepbbarbnemum_rambo, only: ramb

   implicit none

   integer :: NEVT = 1

   integer :: ievt, ierr, prec, particlenr
   real(ki), dimension(8, 4) :: vecs
   real(ki) :: scale2
   real(ki), dimension(0:3) :: amp
   real(ki), dimension(2:3) :: irp
   real(ki) :: t1, t2

   open(unit=logfile,status='unknown',action='write',file='debug.xml')

   open(unit=10,status='old',action='read',file='param.dat',iostat=ierr)
   if(ierr .eq. 0) then
      call parse(10)
      close(unit=10)
   else
      print*, "No file 'param.dat' found. Using defaults"
   end if

   call initgolem()

   call random_seed

   nlo_prefactors=2  ! Do not include any NLO prefactors in order to recognize
                     ! rational numbers for the pole coefficients

   call cpu_time(t1)
   do ievt = 1, NEVT
!      call ramb(5.0E+02_ki**2, vecs)

!      call boost_to_cms(vecs)

!      scale2 = 2.0_ki * dotproduct(vecs(1,:), vecs(2,:))

      vecs(1,:) = (/ 250.0, 0.0, 0.0, 250.0 /)
      vecs(2,:) = (/ 250.0, 0.0, 0.0, -250.0 /)
      vecs(3,:) = (/ 75.44248660279438, -13.21633653678938, -1.177064626981959, -74.26649144952955 /)
      vecs(4,:) = (/ 44.98131949428399, 27.80765763034380, -19.96763010585797, -29.17785167753225 /)
      vecs(5,:) = (/ 92.62223324355639, 44.11745353258120, 79.03574315806995, 19.64382061300453 /)
      vecs(6,:) = (/ 156.4430715983779, -88.78355499018868, -95.32442077834540, 86.63238319307789 /)
      vecs(7,:) = (/ 48.56688843967434, -5.434957739109776, -19.57544900121795, 44.11355452919156 /)
      vecs(8,:) = (/ 81.94400052131301, 35.50973810316283, 57.00882135433334, -46.94541520821219 /)
      
      scale2 = 172.0**2 

! Print generated momenta

      do particlenr = 1, 8
         print*, "Particle ", particlenr, vecs(particlenr,:)
end do


      call print_parameters(scale2)
      call samplitude(vecs, scale2, amp, prec)
      call ir_subtraction(vecs, scale2, irp)
      if(ievt.eq.NEVT) then
         write(*,'(A1,1x,A17,1x,G23.16)') "#", "LO:", amp(0)
         write(*,'(A1,1x,A17,1x,G23.16)') "#", "NLO, finite part:", &
             &    amp(1)/amp(0)
         write(*,'(A1,1x,A17,1x,G23.16)') "#", "NLO, single pole:", &
             &    amp(2)/amp(0)
         write(*,'(A1,1x,A17,1x,G23.16)') "#", "NLO, double pole:", &
             &    amp(3)/amp(0)
      
      
         write(*,'(A1,1x,A17,1x,G23.16)') "#", "IR,  single pole:", &
            & irp(2)/amp(0)
         write(*,'(A1,1x,A17,1x,G23.16)') "#", "IR,  double pole:", &
            & irp(3)/amp(0)
      
      end if
   end do
   call cpu_time(t2)
   write(*,'(A1,1x,A17,1x,F23.3)') "#", "Time/Event [ms]:", &
      & 1.0E+03 * (t2 - t1) / real(NEVT)

   close(logfile)
   call exitgolem()

 contains

subroutine  print_parameters(scale2)
   use p14_ubaru_nmubarepbbarbnemum_config, only: renormalisation, &
        convert_to_cdr, reduction_interoperation, &
        reduction_interoperation_rescue, PSP_check, PSP_rescue
   use p14_ubaru_nmubarepbbarbnemum_model
   implicit none
   real(ki) :: scale2


   write(*,'(A1,1x,A26)') "#", "--------- SETUP ---------"
   write(*,'(A1,1x,A18,I2)') "#", "renormalisation = ", renormalisation

   if(convert_to_cdr) then
      write(*,'(A1,1x,A9,A3)') "#", "scheme = ", "CDR"
   else
      write(*,'(A1,1x,A9,A4)') "#", "scheme = ", "DRED"
   end if

   if(reduction_interoperation.eq.0) then
      write(*,'(A1,1x,A15,A7)') "#", "reduction with ", "SAMURAI"
   else if(reduction_interoperation.eq.1) then
      write(*,'(A1,1x,A15,A7)') "#", "reduction with ", "GOLEM95"
   else if(reduction_interoperation.eq.2) then
      write(*,'(A1,1x,A15,A5)') "#", "reduction with ", "NINJA"
   end if

   if (reduction_interoperation_rescue.ne.reduction_interoperation.and.&
    & PSP_check.and.PSP_rescue) then
      if(reduction_interoperation_rescue.eq.0) then
         write(*,'(A1,1x,A12,A7)') "#", "rescue with ", "SAMURAI"
      else if(reduction_interoperation_rescue.eq.1) then
         write(*,'(A1,1x,A12,A7)') "#", "rescue with ", "GOLEM95"
      else if(reduction_interoperation_rescue.eq.2) then
         write(*,'(A1,1x,A12,A5)') "#", "rescue with ", "NINJA"
      end if
   end if

   write(*,'(A1,1x,A25)') "#", "--- PARAMETER VALUES ---"
   write(*,'(A1,1x,A13)') "#", "Boson masses & widths:"
   write(*,'(A1,1x,A5,G23.16)') "#", "mZ = ", mZ
   write(*,'(A1,1x,A5,G23.16)') "#", "mW = ", mW
   write(*,'(A1,1x,A5,G23.16)') "#", "mH = ", mH
   write(*,'(A1,1x,A5,G23.16)') "#", "wZ = ", wZ
   write(*,'(A1,1x,A5,G23.16)') "#", "wW = ", wW
   write(*,'(A1,1x,A5,G23.16)') "#", "wH = ", wH
   write(*,'(A1,1x,A20)') "#", "Active light quarks:"
   write(*,'(A1,1x,A7,G23.16)') "#", "Nf    =", Nf
   write(*,'(A1,1x,A7,G23.16)') "#", "Nfgen =", Nfgen
   write(*,'(A1,1x,A13)') "#", "Fermion masses & width:"
   write(*,'(A1,1x,A7,G23.16)') "#", "mc   = ", mC
   write(*,'(A1,1x,A7,G23.16)') "#", "mb   = ", mB
   write(*,'(A1,1x,A7,G23.16)') "#", "mbMS = ", mBMS
   write(*,'(A1,1x,A7,G23.16)') "#", "wb   = ", wB
   write(*,'(A1,1x,A7,G23.16)') "#", "mt   = ", mT
   write(*,'(A1,1x,A7,G23.16)') "#", "wt   = ", wT
   write(*,'(A1,1x,A7,G23.16)') "#", "mtau = ", mtau
   write(*,'(A1,1x,A7,G23.16)') "#", "wtau = ", wtau
   write(*,'(A1,1x,A25)') "#", "-------------------------"
   write(*,'(A1,1x,A22)') "#", "Renormalisation scale:"
   write(*,'(A1,1x,A5,G23.16)') "#", "mu = ", sqrt(scale2)
   write(*,'(A1,1x,A25)') "#", "-------------------------"



end subroutine print_parameters

end program test
