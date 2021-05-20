subroutine val(coord_Mo, coord_C, suma, dev, ang, op, cal_eng)
  implicit none
  real(kind=8), intent(in) :: coord_C(576,4), coord_Mo(576,3)
  !
  real(kind=8), intent(out) :: dev, cal_eng
  integer, intent(out) :: suma, ang
  integer :: nc, nmo
  !
  real(kind=8) :: crd(864,3)
  !
  integer :: i,j,m,num(3000),nas,n,suma2(3000), ang_raw
  real(kind=8) :: dx=0.0, dy=0.0, dz=0.0, dist=0.0, sx=0.0, sy=0.0, sz=0.0
  real(kind=8) :: sum_bv_moc(3000), sum_bv_momo(3000), sum_bv_cc(3000), sum_bv_cmo(3000)
  real(kind=8) :: bv=0.0, sum_bv_dev=0.0, op
  real(kind=8) :: pi, c(3), vectores(3000,4)
  !
  real(kind=8) :: alpha=0.0, norma1=0.0, norma2=0.0, normat=0.0, xx=0.0, yy=0.0, zz=0.0, pp=0.0
  !
  pi=4*atan(1.0)
  !
  nc=288
  nmo=576
  !
  !suma: es la sumatoria de los carbonos al cuadrado
  !
  !
  !  open(30,file="struct")          !este file en teoria tiene una estructura completa con un arreglo en particular
  !33 FORMAT(1X,F10.5,1X,F10.5,1X,F10.5)

  do i=1,3000              !inicializando estos arreglos, es decir, metiendo ceros en cada espacio de la memoria.
     sum_bv_moc(i)=0.0
     sum_bv_momo(i)=0.0
     sum_bv_cc(i)=0.0
     sum_bv_cmo(i)=0.0
     num(i)=0
     suma2(i)=0            !esto es para el parametro de orden
     vectores(i,:)=0.0
  enddo
  !
  do i=1,nmo                     !esto estaba comentado
     crd(i,:)=coord_Mo(i,:)    !aca se vuelve turbio el asunto
  enddo                         !miedo, terror, espanto
  !

  m=nmo
  do i=1,nmo
     if (int(coord_C(i,4)).ne.0) then
        m=m+1
        do j=1,3
           crd(m,j)=coord_C(i,j)
        enddo
        !print*, (crd(m,j), j=1,3)
     endif
  enddo

  !print*, m

  !c(1)=0.0036           !C cuadrados from Mo2C - la primera vez con un R2=0.90
  !c(2)=0.1093           !BV from Mo2C          - la primera vez con un R2=0.90
  !c(3)=0.0036           !angles from Mo2C      - la primera vez con un R2=0.90

  c(1)=0.0014*13.605668       !C cuadrados (eV/C^2) |  Mo2C - IB usando R, con un R2=0.957 que se obtuvieron de un sistema 2x2x2
  c(2)=0.6912*13.605668       !BV          (eV/dev) |  Mo2C - IB usando R, con un R2=0.957 que se obtuvieron de un sistema 2x2x2
  c(3)=0.0069*13.605668       !angles      (eV/ang) |  Mo2C - IB usando R, con un R2=0.957 que se obtuvieron de un sistema 2x2x2

  sx=17.0799999237         !lattice x
  sy=25.6300008392         !lattice y
  sz=25.6300008392         !lattice z


  !                               ###########################################################################
  !                               ###########################  carbonos al  #################################
  !                               ###########################   cuadrado    #################################
  !                               ###########################################################################

  do i=1,nmo
     n=0
     do j=nmo+1, nmo+nc
        dx=crd(i,1)-crd(j,1)
        if (dx.gt.sx/2.0) dx=dx-sx
        if (dx.lt.-sx/2.0) dx=dx+sx

        dy=crd(i,2)-crd(j,2)
        if ( dy.gt.(sy/2.0) ) dy=dy-sy
        if ( dy.lt.(-sy/2.0)) dy=dy+sy

        dz=crd(i,3)-crd(j,3)
        if ( dz.gt.(sz/2.0) ) dz=dz-sz
        if ( dz.lt.(-sz/2.0)) dz=dz+sz

        dist=sqrt(dx**2+dy**2+dz**2)

        if (dist.lt.2.4) then
           n=n+1
        endif
     enddo
     num(i)=n                 ! estoy guardando/almacenando cada valor de n en el arreglo num(i), eso es lo que se hace aca.
  enddo

  suma=0
  do i=1,nmo
     suma=suma+num(i)**2     !aca esta haciendo la sumatoria de los carbonos al cuadrado
     suma2(i)=num(i)**2
     !  write(*,*)i,num(i)**2
  enddo

  op=0.0
  op=sqrt(sum(suma2,dim=1)*(1.0/nmo))  !parametro de orden, Root Mean Square (RMS)

  ! write(*,*)suma   !valor final de la sumatoria de los carbonos al cuadrado


  !                               ###########################################################################
  !                               ###########################              ##################################
  !                               ########################### BOND VALENCE ##################################
  !                               ###########################              ##################################
  !                               ###########################################################################


  !Bond valence for Mo-C -- I
  !R is the length of a bond between the two given atoms. In this code R is represented by "dist"
  !The bond valence has the property that its sum around each atom in a compound is equal to the valence (oxidation state) of that atom.
  !RMoC=1.877

  !OPEN(UNIT=3000, FILE="BV_MoC.dat", ACTION="WRITE")
  bv=0.0
  do i=1,nmo                        !este bucle have BV para Mo-C, usando de 1-32 estoy usando los Mo
     !n=0
     bv=0.0
     do j=nmo+1,nmo+nc                    !aca estoy usando los C
        if (i.ne.j) then
           dx=crd(i,1)-crd(j,1)
           if (dx.gt.sx/2.0) dx=dx-sx
           if (dx.lt.-sx/2.0) dx=dx+sx

           dy=crd(i,2)-crd(j,2)
           if (dy.gt.sy/2.0) dy=dy-sy
           if (dy.lt.-sy/2.0 ) dy=dy+sy

           dz=crd(i,3)-crd(j,3)
           if (dz.gt.sz/2.0) dz=dz-sz
           if (dz.lt.-sz/2.0) dz=dz+sz

           dist=sqrt(dx**2+dy**2+dz**2)

           if (dist.lt.2.3) then
              bv=EXP(((1.985-dist)*(1.0/0.37)))    !Mo-C bond-valence parameter
              sum_bv_moc(i)=sum_bv_moc(i)+bv
              !n=n+1
              !write(111,*)i,(j-32)
           endif
        endif
     enddo
     !write(*,*)sum_bv_moc(i),i
     !write(111,*)i,j
  enddo
  !write(*,*)


  !Bond valence for Mo-Mo -- II

  !OPEN(UNIT=1001, FILE="BV_MoMo.dat", ACTION="WRITE")
  !write(1001,*)"    BV MoMo         Mo atom"

  bv=0.0
  do i=1,nmo    !este buble hace BV para Mo-Mo
     !n=0
     bv=0.0
     do j=1,nmo
        if (i.ne.j) then
           dx=crd(i,1)-crd(j,1)
           if (dx.gt.sx/2.0) dx=dx-sx
           if (dx.lt.-sx/2.0) dx=dx+sx

           dy=crd(i,2)-crd(j,2)
           if (dy.gt.sy/2.0) dy=dy-sy
           if (dy.lt.-sy/2.0) dy=dy+sy

           dz=crd(i,3)-crd(j,3)
           if (dz.gt.sz/2.0) dz=dz-sz
           if (dz.lt.-sz/2.0) dz=dz+sz

           dist=sqrt(dx**2+dy**2+dz**2)

           if (dist.lt.3.3) then
              bv=EXP( ((2.615-dist)*(1.0/0.37)) )    !Mo-Mo bond-valence parameter
              sum_bv_momo(i)=sum_bv_momo(i)+bv
              !n=n+1
           endif
        endif
     enddo
     !write(*,*)sum_bv_momo(i),i
  enddo
  !write(*,*)


  !BOND VALENCE FOR C-C -- III
  !RCC=1.540

  !open(unit=1001, file="bv_cc.dat", action="write")
  !write(1001,*)"    BV CC         C atom"

  bv=0.0
  do i=nmo+1,nmo+nc
     !n=0
     bv=0.0
     do j=nmo+1, nmo+nc                      !C atoms
        if (i.ne.j) then
           dx=crd(i,1)-crd(j,1)
           if (dx.gt.sx/2.0) dx=dx-sx
           if (dx.lt.-sx/2.0) dx=dx+sx

           dy=crd(i,2)-crd(j,2)
           if (dy.gt.sy/2.0) dy=dy-sy
           if (dy.lt.-sy/2.0) dy=dy+sy

           dz=crd(i,3)-crd(j,3)
           if (dz.gt.sz/2.0) dz=dz-sz
           if (dz.lt.-sz/2.0) dz=dz+sz

           dist=sqrt(dx**2+dy**2+dz**2)

           if (dist.lt.3.4) then
              bv=EXP( ((1.540-dist)*(1.0/0.37)) )    !C-C bond-valence parameter
              sum_bv_cc(i)=sum_bv_cc(i)+bv
           endif
        endif
     enddo
     !write(*,*) sum_bv_cc(i),i
  enddo
  !write(*,*)

  !BV for C-Mo -- IV

  !open(unit=1001, file="bv_cmo.dat", action="write")
  !write(1001,*)"    BV CMo         Mo atom"

  bv=0.0
  do i=nmo+1,nmo+nc             !atomos de carbono.
     !n=0
     bv=0.0
     do j=1,nmo             !atomos de mo
        if (i.ne.j) then
           dx=crd(i,1)-crd(j,1)
           if (dx.gt.sx/2.0) dx=dx-sx
           if (dx.lt.-sx/2.0) dx=dx+sx

           dy=crd(i,2)-crd(j,2)
           if (dy.gt.sy/2.0) dy=dy-sy
           if (dy.lt.-sy/2.0) dy=dy+sy

           dz=crd(i,3)-crd(j,3)
           if (dz.gt.sz/2.0) dz=dz-sz
           if (dz.lt.-sz/2.0) dz=dz+sz

           dist=sqrt(dx**2+dy**2+dz**2)

           if (dist.lt.2.3) then
              bv=EXP( ((1.985-dist)*(1.0/0.37)) )    !C-Mo bond-valence parameter
              sum_bv_cmo(i)=sum_bv_cmo(i)+bv
           endif
        endif
     enddo
     !write(*,*)sum_bv_cmo(i),i
  enddo
  !write(*,*)

  sum_bv_dev=0.0
  do i=1,nmo
     sum_bv_dev=sum_bv_dev + ((sum_bv_moc(i) + sum_bv_momo(i)-6)**2)             !sum of ALL BV values for Mo.
  enddo                                                                          !cual es el estado de oxidacion en el Mo2C ??

  !write(*,*)sum_bv_dev

  !do i=1,32
  !   write(*,*) sum_bv_moc(i) + sum_bv_momo(i)                                  !aca imprime la valencia de cada Mo individual
  !enddo

  !write(*,*)
  dev=0.0
  dev=sqrt((sum_bv_dev/nmo))  !desviacion promedio de los atomos de Mo.
  !write(*,*)dev
  !write(*,*)sum_bv_dev,dev

  !write (*,*)suma,dev


  !                               ###########################################################################
  !                               ###########################              ##################################
  !                               ###########################  Angles 180  ##################################
  !                               ###########################              ##################################
  !                               ###########################################################################

  nas=0.0                       !el truco para el contador raro, tratar de entender.
  do i=1,nmo                    !almacenando las coordenadas de Molybdeno
     n=0
     xx=0.0
     yy=0.0
     zz=0.0
     do j=nmo+1,nmo+nc               !atomos de Carbono
        dx=crd(i,1)-crd(j,1)
        if (dx.gt.(sx/2.0)) dx=dx-sx
        if (dx.lt.(-sx/2.0)) dx=dx+sx

        dy=crd(i,2)-crd(j,2)
        if (dy.gt.(sy/2.0)) dy=dy-sy
        if (dy.lt.(-sy/2.0)) dy=dy+sy

        dz=crd(i,3)-crd(j,3)
        if (dz.gt.(sz/2.0)) dz=dz-sz
        if (dz.lt.(-sz/2.0)) dz=dz+sz

        dist=sqrt(dx**2+dy**2+dz**2)

        if (dist.lt.2.4) then
           n=n+1
           xx=dx !crd(i,1)-crd(j,1)         !al hacer esta resta de un punto menos otro punto, esto se convierte en un vector (imp. no olvidar).
           yy=dy !crd(i,2)-crd(j,2)         !por eso es que hago la resta de coordenada por coordenada, para obtener cada vector.
           zz=dz !crd(i,3)-crd(j,3)
           !42         FORMAT(F10.1,1X,F10.5,1X,F10.5,1X,F10.5)
           !write(*,42)i,xx,yy,zz
           nas=nas+1
           vectores(nas,:)=(/dble(i), xx, yy, zz/)   !aprender que este es el formato que se usa el / / para encerrar las cantidades deseadas.
           !write(*,*)i, n, vectores(nas,:)
           !write(*,42)vectores(nas,:)
        endif
     enddo
  enddo

  ang=0     !si quito el comentario para 'ang=0' los angulos en el annealing son iguales a cero, asi que no quitar nunca.
  ang_raw=0
  do i=1,3000
     do n=1,3000
        pp=0.0
        alpha=0.0
        norma1=0.0
        norma2=0.0
        normat=0.0
        if (int(vectores(i,1)).eq.int(vectores(n,1))) then
           pp=dot_product((vectores(i,2:4)),(vectores(n,2:4)))
           norma1=sqrt(((vectores(i,2)**2) + (vectores(i,3)**2) + (vectores(i,4)**2)))
           norma2=sqrt(((vectores(n,2)**2) + (vectores(n,3)**2) + (vectores(n,4)**2)))
           normat=(norma1*norma2)
           alpha=(180/pi)*acos((pp/normat))
           !print*,alpha
           if (alpha.gt.171.and.alpha.lt.185) then
              !print*,alpha
              ang_raw=ang_raw+1
           endif
           ang=ang_raw/2
        endif
     enddo
  enddo


  !OPEN(UNIT=3000, FILE="results.dat", ACTION="WRITE")

  cal_eng=0.0
  cal_eng=cal_eng + (c(1)*suma) + (c(2)*dev) + (c(3)*(ang) )
  !cal_eng=cal_eng + (c(1)*suma) + (c(2)*dev) + (c(3)*(ang/0.5) - 4643.7120 )

  !34 FORMAT(i4,1X,F10.5,1X,i6,1X,F10.6)
  !write(*,34)suma,dev,ang,energy
  !write(*,*)
  !write(*,*)suma, dev, ang/2, op
  !write(*,*)

  return
end subroutine val
