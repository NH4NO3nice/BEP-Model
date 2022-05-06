!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!时间差分!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
MODULE module_time_integration
USE module_space_integration           !  声明空间差分模块
    contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!欧拉-后插格式!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!与“三步法”起步的中央差格式组合使用!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!本子程序得到的结果带有边界!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Euler_post(ub, vb, zb, uc, vc, zc, m, f, d, dt, nx, ny, space_integration_type)
    ! ub, vb, zb：当前时刻物理量
    ! uc, vc, zc：下一时刻物理量
    ! E, G, H：物理量倾向
    ! m, f：地图放大系数、地转参数
    ! d, dt：空间步长，时间步长
    ! nx, ny：网格格点
    ! gt：重力加速度
    ! space_integration_type：空间差分格式
        IMPLICIT NONE
        INTEGER :: i, j
        INTEGER :: nx, ny, space_integration_type, dt
        REAL, PARAMETER :: gt = 9.8
        REAL :: d
        REAL :: ub(nx, ny), vb(nx, ny), zb(nx, ny)
        REAL :: uc(nx, ny), vc(nx, ny), zc(nx, ny)
        REAL :: E(nx, ny), G(nx, ny), H(nx, ny)
        REAL :: m(nx, ny), f(nx, ny)
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!计算本时刻预报量倾向!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (space_integration_type == 1) then
            call advection_equation_2(ub, vb, zb, m, f, E, G, H, d, nx, ny)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!下一时刻初报结果，即带星号的uvz!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i = 2, nx - 1
            do j = 2, ny - 1
                uc(i, j) = ub(i, j) + dt * E(i, j)
                vc(i, j) = vb(i, j) + dt * G(i, j)
                zc(i, j) = zb(i, j) + dt * H(i, j)
            enddo
        enddo
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!计算带星号的预报量倾向!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !call fixed_boundary(ub, vb, zb, uc, vc, zc, nx, ny)  !  为带星号的预报量添加边界条件
        if (space_integration_type == 1) then
            call advection_equation_2(uc, vc, zc, m, f, E, G, H, d, nx, ny)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!补充预报下一时刻!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i = 2, nx - 1
            do j = 2, ny - 1
                uc(i, j) = ub(i, j) + dt * E(i, j)
                vc(i, j) = vb(i, j) + dt * G(i, j)
                zc(i, j) = zb(i, j) + dt * H(i, j)
            enddo
        enddo
        return
    END SUBROUTINE Euler_post    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!中央差“三步法”起报!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!在欧拉-后差预报的基础上“三步法”起报，向前报两步!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!减小中央差数值解的影响!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!本程序计算包括边界!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE three_step_start(ua, va, za, ub, vb, zb, uc, vc, zc, m, f, d, dt, nx, ny, space_integration_type)
    ! ua, va, za：前一时刻物理量
    ! ub, vb, zb：当前时刻物理量
    ! uc, vc, zc：下一时刻物理量
    ! E, G, H：物理量倾向
    ! m, f：地图放大系数、地转参数
    ! d, dt：空间步长，时间步长
    ! nx, ny：网格格点
    ! gt：重力加速度
    ! space_integration_type：空间差分格式
        IMPLICIT NONE
        INTEGER :: i, j
        INTEGER :: nx, ny, space_integration_type, dt
        REAL, PARAMETER :: gt = 9.8
        REAL :: d
        REAL :: ua(nx, ny), va(nx, ny), za(nx, ny)  ! 用于储存第二步计算的结果，为后面中央差分创造条件
        REAL :: ub(nx, ny), vb(nx, ny), zb(nx, ny)
        REAL :: uc(nx, ny), vc(nx, ny), zc(nx, ny)
        REAL :: E(nx, ny), G(nx, ny), H(nx, ny)
        REAL :: m(nx, ny), f(nx, ny)
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!第一步：向前半步差分!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (space_integration_type == 1) then
            call advection_equation_2(ub, vb, zb, m, f, E, G, H, d, nx, ny)
         endif   !  计算本时刻预报量倾向
     
        do i = 2, nx - 1
            do j = 2, ny - 1
                uc(i, j) = ub(i, j) + dt * E(i, j) / 2.0
                vc(i, j) = vb(i, j) + dt * G(i, j) / 2.0
                zc(i, j) = zb(i, j) + dt * H(i, j) / 2.0
            enddo
        enddo
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!第二步：半步中央差分!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !call fixed_boundary(ub, vb, zb, uc, vc, zc, nx, ny)  !  添加边界条件
        if (space_integration_type == 1) then
            call advection_equation_2(uc, vc, zc, m, f, E, G, H, d, nx, ny)
        endif   !  计算本时刻预报量倾向
     
        do i = 2, nx - 1
            do j = 2, ny - 1
                uc(i, j) = ub(i, j) + dt * E(i, j)
                vc(i, j) = vb(i, j) + dt * G(i, j)
                zc(i, j) = zb(i, j) + dt * H(i, j)
            enddo
        enddo
    
        call transmit(uc, vc, zc, ua, va, za, nx, ny)  !  将第二步的预报结果储存在ua,va,za当中
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!第三步：中央差分!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  第二步并没有改变边界值，不需要重新传递
        if (space_integration_type == 1) then
            call advection_equation_2(uc, vc, zc, m, f, E, G, H, d, nx, ny)
        endif   !  计算本时刻预报量倾向
     
        do i = 2, nx - 1
            do j = 2, ny - 1
                uc(i, j) = ub(i, j) + 2.0 * dt * E(i, j)
                vc(i, j) = vb(i, j) + 2.0 * dt * G(i, j)
                zc(i, j) = zb(i, j) + 2.0 * dt * H(i, j)
            enddo
        enddo
        return
        END SUBROUTINE three_step_start
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!中央差!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!在欧拉-后差、“三步法”起报的基础上!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!本程序只计算内网格!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE central_difference(ua, va, za, ub, vb, zb, uc, vc, zc, m, f, d, dt, nx, ny, space_integration_type)
    ! ua, va, za：前一时刻物理量
    ! ub, vb, zb：当前时刻物理量
    ! uc, vc, zc：下一时刻物理量
    ! E, G, H：物理量倾向
    ! m, f：地图放大系数、地转参数
    ! d, dt：空间步长，时间步长
    ! nx, ny：网格格点
    ! gt：重力加速度
    ! space_integration_type：空间差分格式
        IMPLICIT NONE
        INTEGER :: i, j
        INTEGER :: nx, ny, space_integration_type, dt
        REAL, PARAMETER :: gt = 9.8
        REAL :: d
        REAL :: ua(nx, ny), va(nx, ny), za(nx, ny)
        REAL :: ub(nx, ny), vb(nx, ny), zb(nx, ny)
        REAL :: uc(nx, ny), vc(nx, ny), zc(nx, ny)
        REAL :: E(nx, ny), G(nx, ny), H(nx, ny)
        REAL :: m(nx, ny), f(nx, ny)
    
        if (space_integration_type == 1) then
            call advection_equation_2(ub, vb, zb, m, f, E, G, H, d, nx, ny)
        endif   !  计算本时刻预报量倾向
     
        do i = 2, nx - 1
            do j = 2, ny - 1
                uc(i, j) = ua(i, j) + 2.0 * dt * E(i, j)
                vc(i, j) = va(i, j) + 2.0 * dt * G(i, j)
                zc(i, j) = za(i, j) + 2.0 * dt * H(i, j)
            enddo
        enddo
    END SUBROUTINE central_difference
END MODULE
