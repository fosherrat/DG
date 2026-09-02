subroutine update_artificial_viscosity()
    use m_parameter
    use m_input
    use m_point_dg
    use m_flow

    implicit none

    integer :: elem, face, k, q, vl, vr, npoint
    double precision :: energy, energy_p, ratio, s0, eps0, eps_face
    double precision :: xi, rho, vel, pressure, speed, max_speed
    double precision :: u(nvar)
    double precision, external :: poly_value

    if (dim.ne.1) then
        write(*,*) 'ERROR: artificial viscosity currently supports only dim = 1'
        stop
    endif

    if (.not.allocated(shock_sensor)) allocate(shock_sensor(nvol))
    if (.not.allocated(av_elem)) allocate(av_elem(nvol))
    if (.not.allocated(av_face)) allocate(av_face(nsur_tot))
    if (.not.allocated(av_end)) allocate(av_end(2,nvol))

    shock_sensor = -huge(1.d0)
    av_elem = 0.d0
    av_face = 0.d0
    av_end = 0.d0

    if (shock_capture.eq.0) return

    s0 = -4.d0*log10(dble(order))
    npoint = max(ndof+1,3)

    do elem=1,nvol
        energy = 0.d0
        do k=1,ndof
            energy = energy + cons(1,k,elem)**2*mass(k,k)
        enddo
        energy_p = cons(1,ndof,elem)**2*mass(ndof,ndof)
        ratio = energy_p/max(energy,tiny(1.d0))
        shock_sensor(elem) = log10(max(ratio,1.d-30))

        max_speed = 0.d0
        do q=1,npoint
            xi = -1.d0 + 2.d0*dble(q-1)/dble(npoint-1)
            u = 0.d0
            do k=1,ndof
                u(:) = u(:) + cons(:,k,elem)*poly_value(legendre(:,k),ndof,xi)
            enddo

            rho = u(1)
            if (rho.le.0.d0) then
                write(*,*) 'ERROR: non-positive density in artificial viscosity sensor'
                stop
            endif
            vel = u(2)/rho
            pressure = (gamma-1.d0)*(u(3)-0.5d0*rho*vel*vel)
            if (pressure.le.0.d0) then
                write(*,*) 'ERROR: non-positive pressure in artificial viscosity sensor'
                stop
            endif
            speed = abs(vel) + sqrt(gamma*pressure/rho)
            max_speed = max(max_speed,speed)
        enddo

        eps0 = av_c*vol(elem)*max_speed/dble(order)
        if (shock_sensor(elem).le.s0-av_kappa) then
            av_elem(elem) = 0.d0
        elseif (shock_sensor(elem).ge.s0+av_kappa) then
            av_elem(elem) = eps0
        else
            av_elem(elem) = 0.5d0*eps0*(1.d0 &
                + sin(pi*(shock_sensor(elem)-s0)/(2.d0*av_kappa)))
        endif
    enddo

    ! C0 regularization in 1D: use the maximum neighboring element viscosity
    ! at each mesh vertex, then interpolate linearly inside each element.
    do face=1,nsur_tot
        vl = con_sur_vol(1,face)
        vr = con_sur_vol(2,face)

        eps_face = av_elem(vl)
        if (vr.gt.0) eps_face = max(eps_face,av_elem(vr))
        av_face(face) = eps_face

        xi = (sur_cen(1,face)-vol_cen(1,vl))/jcb(1,1,vl)
        if (xi.lt.0.d0) then
            av_end(1,vl) = max(av_end(1,vl),eps_face)
        else
            av_end(2,vl) = max(av_end(2,vl),eps_face)
        endif

        if (vr.gt.0) then
            xi = (sur_cen(1,face)-vol_cen(1,vr))/jcb(1,1,vr)
            if (xi.lt.0.d0) then
                av_end(1,vr) = max(av_end(1,vr),eps_face)
            else
                av_end(2,vr) = max(av_end(2,vr),eps_face)
            endif
        endif
    enddo

end subroutine update_artificial_viscosity
