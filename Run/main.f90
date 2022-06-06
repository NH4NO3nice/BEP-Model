!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!正压原始方程模式!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM main
! 纬向格点数nx
! 经向格点数ny  
! 区域中心经纬度clon, clat
! 网格距d/米
! 积分步长dt /秒
! 地图放大系数m，地转参数f
! n-1，n，n+1时间层的位势高度场、u风场、v风场za(:, :), ua(:, :), va(:, :), zb(:, :), ub(:, :), vb(:, :), zc(:, :), uc(:, :), vc(:, :)
! 初边条件initial, boundary
! 时间差分格式time_integration_type
! 空间差分格式space_integration_type
! 总积分时长all_time/秒
! 时空平滑系数S_time, S_space
! 输出数据时间间隔output_interval/秒 
! 存放初始时间start_time
! 是否做正逆平滑 Plus_minus_smooth
! 定义初始数据文件地址，文件地址及预报结果文件名 input_path, output_path, output_filename

USE module_initialization                   !  声明初始化模块
USE module_boundary                        !  声明边界条件模块
USE module_m_f                                   !  声明放大系数和地转参数计算模块
USE module_time_integration             !  声明时间差分模块
USE module_space_smooth                 !  声明空间平滑模块
USE module_time_smooth                   !  声明时间平滑模块

    IMPLICIT NONE
    
    INTEGER :: nx, ny  !  定义区域网格点
    INTEGER :: i, j, p, q  !  循环变量
    INTEGER :: initial, boundary, time_integration_type, space_integration_type
    INTEGER :: Plus_minus_smooth  ! 是否做正逆平滑
    INTEGER :: all_time, dt  !  定义总积分时长, 积分步长
    INTEGER :: output_interval  !  输出数据时间间隔/秒
    INTEGER :: int, remainder  !  all_time与12h相除的除数和余数
    REAL :: d  !  定义网格距
    REAL :: clon, clat  !  定义区域中心经纬度
    REAL :: S_time, S_space
    REAL, ALLOCATABLE :: m(:, :), f(:, :)  !  定义地图放大系数m和地转参数f的数组，可变数组
    REAL, ALLOCATABLE :: za(:, :), ua(:, :), va(:, :), &
                                               zb(:, :), ub(:, :), vb(:, :), &
                                               zc(:, :), uc(:, :), vc(:, :)  !  定义n-1，n，n+1时间层的位势高度场、u风场、v风场
    CHARACTER(50) :: input_path, output_path, output_filename  !  定义初始数据文件地址，文件地址及预报结果文件名
    CHARACTER :: start_time  !  存放初始时间
    
    write(*, '(100A)'), '################################################################'
    write(*, '(100A)'), '     Welcome to                      '
    write(*, '(100A)'), '     USE THE BAROTROPIC PRIMITIVE EQUSTION MODEL '
    write(*, '(100A)'), '     By Nanjing University of Information and Technology      '
    write(*, '(100A)'), '     School of Atomospheric Physics  '
    write(*, '(100A)'), '     Shanchuan Xiao    '
    write(*, '(100A)'), '################################################################'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 从namelist文件中读取控制变量!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    namelist / control / start_time, all_time, output_interval, nx, ny, d, dt, clon, clat, input_path, output_path, output_filename
    namelist / initial_boundary / initial, boundary
    namelist / integration / time_integration_type, space_integration_type
    namelist / smooth / Plus_minus_smooth, S_time, S_space
    
    open(11, file = 'C:\Users\NH4NO3nice\Desktop\PE_module\namelist')
    read(11, nml = control)
    open(22, file = 'C:\Users\NH4NO3nice\Desktop\PE_module\namelist')
    read(22, nml = initial_boundary)
    open(33, file = 'C:\Users\NH4NO3nice\Desktop\PE_module\namelist')
    read(33, nml = integration)
    open(44, file = 'C:\Users\NH4NO3nice\Desktop\PE_module\namelist')
    read(44, nml = smooth)
    close(11); close(22); close(33); close(44)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 读入位势高度场数据!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(za(nx, ny))  ! 为位势高度场数组分配内存
    open(12, file = trim(input_path) // 'za.grd', form = 'binary') 
    read(12), ((za(i, j), i = 1, nx), j = 1, ny)
    close(12)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 计算区域网格上的地图放大系数m和地转参数f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(m(nx, ny), f(nx, ny))  !  为可变数组分配内存
    !m = 0; f = 0
    call m_f(m, f, d, clat, nx, ny)  !  调用m_f外部子例行程序
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输出m，f数据!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(15, file = trim(output_path) // 'm.txt')
    open(16, file = trim(output_path) // 'f.txt')
    write(15, '(20f10.5)') m
    write(16, '(20f10.5)') f
    close(15); close(16)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 初始化!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(ua(nx, ny), va(nx, ny))
    !ua = 0; va = 0
    if (initial == 0) then  !  无初始化条件
        open(13, file = trim(input_path) // 'ua.grd', form = 'binary')
        open(14, file = trim(input_path) // 'va.grd', form = 'binary')
        read(13), ((ua(i, j), i = 1, nx), j = 1, ny)
        read(14), ((va(i, j), i = 1, nx), j = 1, ny)
        close(13);close(14)
        
    else if (initial == 1) then !  中央差地转风初始化
        call geostrophic_wind_centre(ua, va, za, m, f, nx, ny, d)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输出地转风!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(17, file = trim(output_path) // 'geo_u.txt')
        open(18, file = trim(output_path) // 'geo_v.txt')
        write(17, '(20f10.5)'), ua
        write(18, '(20f10.5)'), va
        close(17); close(18)
    else if (initial == 2) then !  前差地转风初始化
        call geostrophic_wind_forward(ua, va, za, m, f, nx, ny, d)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输出地转风!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(17, file = trim(output_path) // 'geo_u.txt')
        open(18, file = trim(output_path) // 'geo_v.txt')
        write(17, '(20f10.5)'), ua
        write(18, '(20f10.5)'), va
        close(17); close(18)
    else if (initial == 3) then !  后差地转风初始化
        call geostrophic_wind_backward(ua, va, za, m, f, nx, ny, d)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输出地转风!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(17, file = trim(output_path) // 'geo_u.txt')
        open(18, file = trim(output_path) // 'geo_v.txt')
        write(17, '(20f10.5)'), ua
        write(18, '(20f10.5)'), va
        close(17); close(18)
    endif
    
    !!!!!!!!!!!!!!!!!!!!!!!将u、v，z初值存入grd中，与最后的预测值放在同一个文件中，方便绘图!!!!!!!!!!!!!!!!!!!!!!!
    open(19, file = trim(output_path) // trim(output_filename), form = 'binary')
    write(19), ((ua(i, j), i= 1, nx), j = 1, ny)
    write(19), ((va(i, j), i= 1, nx), j = 1, ny)
    write(19), ((za(i, j), i= 1, nx), j = 1, ny)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!数值计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(ub(nx, ny), vb(nx, ny), zb(nx, ny))  !  为本时刻物理量分配内存
    allocate(uc(nx, ny), vc(nx, ny), zc(nx, ny))   !  为下一时刻物理量分配内存
    !ub = 0; vb = 0; zb = 0
    !uc = 0; vc = 0; zc = 0
    call transmit(ub, vb, zb, ua, va, za, nx, ny)   !  将初始值传递给当前时刻物理量
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!添加边界条件!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (boundary == 1) then   ! 固定边界
        call fixed_boundary(ub, vb, zb, uc, vc, zc, nx, ny)  
    endif
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!时间差分!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (time_integration_type == 1) then  !  欧拉-后差 + “三步法”起步中央差
        int = all_time / (12 * 3600); remainder = mod(all_time, (12 * 3600))
        if (remainder > 0) then
            int = int + 1
            remainder = remainder / dt  !  将余数(s)转换成积分次数
        endif
        do i = 1, int  !  每12h循环一次 
            print*,i
            do j = 1, 3600 / dt  !  前1h用欧拉-后差格式
                call Euler_post(ub, vb, zb, uc, vc, zc, m, f, d, dt, nx, ny, space_integration_type)  ! 时间差分
                call transmit(ub, vb, zb, uc, vc, zc, nx, ny)  !  将uc, vc, zc的值传递给ub, vb, zb，方便下一步的时间积分
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!判断是否积分终止!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (i == int .and. j == remainder) then
                    print*,j,'stop'
                    write(19), ((ub(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((vb(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((zb(p, q), p= 1, nx), q = 1, ny)
                    close(19)
                    stop  !  终止程序 
                endif
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输出预报结果!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (mod((i - 1) * 12 * 3600 + j * dt, output_interval) == 0) then
                    write(19), ((ub(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((vb(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((zb(p, q), p= 1, nx), q = 1, ny)
                endif
            enddo
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 边界平滑 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call nine_smooth_space_out(ub, nx, ny, S_space)
            call nine_smooth_space_out(vb, nx, ny, S_space)
            call nine_smooth_space_out(zb, nx, ny, S_space)
            
            if (Plus_minus_smooth == 0) then  ! 不做正逆平滑，只做正平滑
                continue  ! 程序继续向下执行；相当于python中的pass
                
            else if (Plus_minus_smooth == 1) then  ! 做正逆平滑
                call nine_smooth_space_out(ub, nx, ny, -1 * S_space)  ! 逆平滑
                call nine_smooth_space_out(vb, nx, ny, -1 * S_space)
                call nine_smooth_space_out(zb, nx, ny, -1 * S_space)
            endif
                
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!“三步法”起报!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call three_step_start(ua, va, za, ub, vb, zb, uc, vc, zc, m, f, d, dt, nx, ny, space_integration_type)  !  “三步法”起步中央差，向前报2步
            call transmit(ub, vb, zb, uc, vc, zc, nx, ny)  !  将uc, vc, zc的值传递给ub, vb, zb，方便下一步的时间积分
            
            do j = 3600 / dt  + 3, 12 * 3600 / dt  !  中央差分12h
                call central_difference(ua, va, za, ub, vb, zb, uc, vc, zc, m, f, d, dt, nx, ny, space_integration_type)  !  中央差
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!时间三点平滑!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (mod(j * dt, 6 * 3600) == 0) then  !  每6h时间平滑一次
                    ! 本步只平滑了n时刻变量，对n+1时刻并未处理，即在后面的中央差用到的两个时次中，一个平滑一个未平滑
                    call three_smooth_time(ua, va, za, ub, vb, zb, uc, vc, zc, nx, ny, S_time)  !  正时间平滑
                    call three_smooth_time(ua, va, za, ub, vb, zb, uc, vc, zc, nx, ny, -1 * S_time)  !  逆时间平滑
                endif
                
                !  不用为中央差结果赋边界值，前面已经有边界了
                call transmit(ua, va, za,ub, vb, zb, nx, ny)  !  将ub, vb, zb的值传递给ua, va, za，方便下一步的时间积分
                call transmit(ub, vb, zb, uc, vc, zc, nx, ny)  !  将uc, vc, zc的值传递给ub, vb, zb，方便下一步的时间积分
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!空间平滑!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (j * dt == 12 * 3600) then  ! 每12h内部网格平滑一次
                    call five_smooth_space_in(ub, nx, ny, S_space)
                    call five_smooth_space_in(vb, nx, ny, S_space)
                    call five_smooth_space_in(zb, nx, ny, S_space)  !  只平滑了当前时间层，n-1时间层并未平滑！！
                    
                    if (Plus_minus_smooth == 0) then  ! 不做正逆平滑，只做正平滑
                        continue  ! 程序继续向下执行；相当于python中的pass
                        
                    else if (Plus_minus_smooth == 1) then
                        call five_smooth_space_in(ub, nx, ny, -1 * s_space)
                        call five_smooth_space_in(vb, nx, ny, -1 * s_space)
                        call five_smooth_space_in(zb, nx, ny, -1 * s_space)  ! 逆平滑
                    endif
                else 
                    if (mod(j * dt, 3600) == 0) then  !  每1h边界平滑一次
                        call nine_smooth_space_out(ub, nx, ny, S_space)
                        call nine_smooth_space_out(vb, nx, ny, S_space)
                        call nine_smooth_space_out(zb, nx, ny, S_space)
                        
                        if (Plus_minus_smooth == 0) then  ! 不做正逆平滑，只做正平滑
                            continue  ! 程序继续向下执行；相当于python中的pass
                            
                        else if (Plus_minus_smooth == 1) then  ! 做正逆平滑
                            call nine_smooth_space_out(ub, nx, ny, -1 * S_space)  ! 逆平滑
                            call nine_smooth_space_out(vb, nx, ny, -1 * S_space)
                            call nine_smooth_space_out(zb, nx, ny, -1 * S_space)
                        endif
                    endif
                endif
                
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!判断是否积分终止!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (i == int .and. j == remainder) then
                    print*, j, 'stop'
                    write(19), ((ub(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((vb(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((zb(p, q), p= 1, nx), q = 1, ny)
                    close(19)
                    stop  !  终止程序 
                endif
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输出预报结果!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (mod((i - 1) * 12 * 3600 + j * dt, output_interval) == 0) then
                    print*, j
                    write(19), ((ub(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((vb(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((zb(p, q), p= 1, nx), q = 1, ny)
                endif
            enddo
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!关闭输出文件!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        close(19)
    endif
END
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!数组值传递子程序!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!便于循环利用已定义的三个时间层的物理量!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE transmit(ua, va, za, ub, vb, zb, nx, ny)
! 将ub, vb, zb的值分别传递给ua, va, za
    IMPLICIT NONE
    INTEGER :: i, j
    INTEGER :: nx, ny
    REAL :: ua(nx, ny), va(nx, ny), za(nx, ny)
    REAL :: ub(nx, ny), vb(nx, ny), zb(nx, ny)
    
    do i = 1, nx
        do j = 1, ny
            ua(i, j) = ub(i, j)
            va(i, j) = vb(i, j)
            za(i, j) = zb(i, j)
        enddo
    enddo
    return
END SUBROUTINE transmit    
    
    

    
    
    