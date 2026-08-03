program dumpr4
  implicit none
  character(len=2) :: forst
  character(len=10) :: voleq
  integer :: spec, errflag, fnum, i
  integer :: specs(30)
  data specs /15,17,19,64,65,66,93,96,101,106,108,113,122,133,142,202, &
              242,299,313,321,322,475,746,748,749,814,998,119,73,263/
  do fnum=1,19
    write(forst,'(I2.2)') fnum
    do i=1,30
      spec=specs(i)
      voleq='          '
      errflag=0
      call R4_EQN(forst, spec, voleq, errflag)
      write(*,'(A,I3,A,I4,A,A)') 'F',fnum,' SP',spec,' ',voleq
    end do
  end do
end program
