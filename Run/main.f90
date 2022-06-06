!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ѹԭʼ����ģʽ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM main
! γ������nx
! ��������ny  
! �������ľ�γ��clon, clat
! �����d/��
! ���ֲ���dt /��
! ��ͼ�Ŵ�ϵ��m����ת����f
! n-1��n��n+1ʱ����λ�Ƹ߶ȳ���u�糡��v�糡za(:, :), ua(:, :), va(:, :), zb(:, :), ub(:, :), vb(:, :), zc(:, :), uc(:, :), vc(:, :)
! ��������initial, boundary
! ʱ���ָ�ʽtime_integration_type
! �ռ��ָ�ʽspace_integration_type
! �ܻ���ʱ��all_time/��
! ʱ��ƽ��ϵ��S_time, S_space
! �������ʱ����output_interval/�� 
! ��ų�ʼʱ��start_time
! �Ƿ�������ƽ�� Plus_minus_smooth
! �����ʼ�����ļ���ַ���ļ���ַ��Ԥ������ļ��� input_path, output_path, output_filename

USE module_initialization                   !  ������ʼ��ģ��
USE module_boundary                        !  �����߽�����ģ��
USE module_m_f                                   !  �����Ŵ�ϵ���͵�ת��������ģ��
USE module_time_integration             !  ����ʱ����ģ��
USE module_space_smooth                 !  �����ռ�ƽ��ģ��
USE module_time_smooth                   !  ����ʱ��ƽ��ģ��

    IMPLICIT NONE
    
    INTEGER :: nx, ny  !  �������������
    INTEGER :: i, j, p, q  !  ѭ������
    INTEGER :: initial, boundary, time_integration_type, space_integration_type
    INTEGER :: Plus_minus_smooth  ! �Ƿ�������ƽ��
    INTEGER :: all_time, dt  !  �����ܻ���ʱ��, ���ֲ���
    INTEGER :: output_interval  !  �������ʱ����/��
    INTEGER :: int, remainder  !  all_time��12h����ĳ���������
    REAL :: d  !  ���������
    REAL :: clon, clat  !  �����������ľ�γ��
    REAL :: S_time, S_space
    REAL, ALLOCATABLE :: m(:, :), f(:, :)  !  �����ͼ�Ŵ�ϵ��m�͵�ת����f�����飬�ɱ�����
    REAL, ALLOCATABLE :: za(:, :), ua(:, :), va(:, :), &
                                               zb(:, :), ub(:, :), vb(:, :), &
                                               zc(:, :), uc(:, :), vc(:, :)  !  ����n-1��n��n+1ʱ����λ�Ƹ߶ȳ���u�糡��v�糡
    CHARACTER(50) :: input_path, output_path, output_filename  !  �����ʼ�����ļ���ַ���ļ���ַ��Ԥ������ļ���
    CHARACTER :: start_time  !  ��ų�ʼʱ��
    
    write(*, '(100A)'), '################################################################'
    write(*, '(100A)'), '     Welcome to                      '
    write(*, '(100A)'), '     USE THE BAROTROPIC PRIMITIVE EQUSTION MODEL '
    write(*, '(100A)'), '     By Nanjing University of Information and Technology      '
    write(*, '(100A)'), '     School of Atomospheric Physics  '
    write(*, '(100A)'), '     Shanchuan Xiao    '
    write(*, '(100A)'), '################################################################'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ��namelist�ļ��ж�ȡ���Ʊ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ����λ�Ƹ߶ȳ�����!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(za(nx, ny))  ! Ϊλ�Ƹ߶ȳ���������ڴ�
    open(12, file = trim(input_path) // 'za.grd', form = 'binary') 
    read(12), ((za(i, j), i = 1, nx), j = 1, ny)
    close(12)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! �������������ϵĵ�ͼ�Ŵ�ϵ��m�͵�ת����f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(m(nx, ny), f(nx, ny))  !  Ϊ�ɱ���������ڴ�
    !m = 0; f = 0
    call m_f(m, f, d, clat, nx, ny)  !  ����m_f�ⲿ�����г���
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���m��f����!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(15, file = trim(output_path) // 'm.txt')
    open(16, file = trim(output_path) // 'f.txt')
    write(15, '(20f10.5)') m
    write(16, '(20f10.5)') f
    close(15); close(16)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ��ʼ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(ua(nx, ny), va(nx, ny))
    !ua = 0; va = 0
    if (initial == 0) then  !  �޳�ʼ������
        open(13, file = trim(input_path) // 'ua.grd', form = 'binary')
        open(14, file = trim(input_path) // 'va.grd', form = 'binary')
        read(13), ((ua(i, j), i = 1, nx), j = 1, ny)
        read(14), ((va(i, j), i = 1, nx), j = 1, ny)
        close(13);close(14)
        
    else if (initial == 1) then !  ������ת���ʼ��
        call geostrophic_wind_centre(ua, va, za, m, f, nx, ny, d)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����ת��!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(17, file = trim(output_path) // 'geo_u.txt')
        open(18, file = trim(output_path) // 'geo_v.txt')
        write(17, '(20f10.5)'), ua
        write(18, '(20f10.5)'), va
        close(17); close(18)
    else if (initial == 2) then !  ǰ���ת���ʼ��
        call geostrophic_wind_forward(ua, va, za, m, f, nx, ny, d)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����ת��!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(17, file = trim(output_path) // 'geo_u.txt')
        open(18, file = trim(output_path) // 'geo_v.txt')
        write(17, '(20f10.5)'), ua
        write(18, '(20f10.5)'), va
        close(17); close(18)
    else if (initial == 3) then !  ����ת���ʼ��
        call geostrophic_wind_backward(ua, va, za, m, f, nx, ny, d)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����ת��!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(17, file = trim(output_path) // 'geo_u.txt')
        open(18, file = trim(output_path) // 'geo_v.txt')
        write(17, '(20f10.5)'), ua
        write(18, '(20f10.5)'), va
        close(17); close(18)
    endif
    
    !!!!!!!!!!!!!!!!!!!!!!!��u��v��z��ֵ����grd�У�������Ԥ��ֵ����ͬһ���ļ��У������ͼ!!!!!!!!!!!!!!!!!!!!!!!
    open(19, file = trim(output_path) // trim(output_filename), form = 'binary')
    write(19), ((ua(i, j), i= 1, nx), j = 1, ny)
    write(19), ((va(i, j), i= 1, nx), j = 1, ny)
    write(19), ((za(i, j), i= 1, nx), j = 1, ny)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ֵ����!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(ub(nx, ny), vb(nx, ny), zb(nx, ny))  !  Ϊ��ʱ�������������ڴ�
    allocate(uc(nx, ny), vc(nx, ny), zc(nx, ny))   !  Ϊ��һʱ�������������ڴ�
    !ub = 0; vb = 0; zb = 0
    !uc = 0; vc = 0; zc = 0
    call transmit(ub, vb, zb, ua, va, za, nx, ny)   !  ����ʼֵ���ݸ���ǰʱ��������
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ӱ߽�����!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (boundary == 1) then   ! �̶��߽�
        call fixed_boundary(ub, vb, zb, uc, vc, zc, nx, ny)  
    endif
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ʱ����!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (time_integration_type == 1) then  !  ŷ��-��� + �����������������
        int = all_time / (12 * 3600); remainder = mod(all_time, (12 * 3600))
        if (remainder > 0) then
            int = int + 1
            remainder = remainder / dt  !  ������(s)ת���ɻ��ִ���
        endif
        do i = 1, int  !  ÿ12hѭ��һ�� 
            print*,i
            do j = 1, 3600 / dt  !  ǰ1h��ŷ��-����ʽ
                call Euler_post(ub, vb, zb, uc, vc, zc, m, f, d, dt, nx, ny, space_integration_type)  ! ʱ����
                call transmit(ub, vb, zb, uc, vc, zc, nx, ny)  !  ��uc, vc, zc��ֵ���ݸ�ub, vb, zb��������һ����ʱ�����
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�ж��Ƿ������ֹ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (i == int .and. j == remainder) then
                    print*,j,'stop'
                    write(19), ((ub(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((vb(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((zb(p, q), p= 1, nx), q = 1, ny)
                    close(19)
                    stop  !  ��ֹ���� 
                endif
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���Ԥ�����!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (mod((i - 1) * 12 * 3600 + j * dt, output_interval) == 0) then
                    write(19), ((ub(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((vb(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((zb(p, q), p= 1, nx), q = 1, ny)
                endif
            enddo
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! �߽�ƽ�� !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call nine_smooth_space_out(ub, nx, ny, S_space)
            call nine_smooth_space_out(vb, nx, ny, S_space)
            call nine_smooth_space_out(zb, nx, ny, S_space)
            
            if (Plus_minus_smooth == 0) then  ! ��������ƽ����ֻ����ƽ��
                continue  ! �����������ִ�У��൱��python�е�pass
                
            else if (Plus_minus_smooth == 1) then  ! ������ƽ��
                call nine_smooth_space_out(ub, nx, ny, -1 * S_space)  ! ��ƽ��
                call nine_smooth_space_out(vb, nx, ny, -1 * S_space)
                call nine_smooth_space_out(zb, nx, ny, -1 * S_space)
            endif
                
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!������������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call three_step_start(ua, va, za, ub, vb, zb, uc, vc, zc, m, f, d, dt, nx, ny, space_integration_type)  !  ������������������ǰ��2��
            call transmit(ub, vb, zb, uc, vc, zc, nx, ny)  !  ��uc, vc, zc��ֵ���ݸ�ub, vb, zb��������һ����ʱ�����
            
            do j = 3600 / dt  + 3, 12 * 3600 / dt  !  ������12h
                call central_difference(ua, va, za, ub, vb, zb, uc, vc, zc, m, f, d, dt, nx, ny, space_integration_type)  !  �����
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ʱ������ƽ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (mod(j * dt, 6 * 3600) == 0) then  !  ÿ6hʱ��ƽ��һ��
                    ! ����ֻƽ����nʱ�̱�������n+1ʱ�̲�δ�������ں����������õ�������ʱ���У�һ��ƽ��һ��δƽ��
                    call three_smooth_time(ua, va, za, ub, vb, zb, uc, vc, zc, nx, ny, S_time)  !  ��ʱ��ƽ��
                    call three_smooth_time(ua, va, za, ub, vb, zb, uc, vc, zc, nx, ny, -1 * S_time)  !  ��ʱ��ƽ��
                endif
                
                !  ����Ϊ����������߽�ֵ��ǰ���Ѿ��б߽���
                call transmit(ua, va, za,ub, vb, zb, nx, ny)  !  ��ub, vb, zb��ֵ���ݸ�ua, va, za��������һ����ʱ�����
                call transmit(ub, vb, zb, uc, vc, zc, nx, ny)  !  ��uc, vc, zc��ֵ���ݸ�ub, vb, zb��������һ����ʱ�����
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�ռ�ƽ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (j * dt == 12 * 3600) then  ! ÿ12h�ڲ�����ƽ��һ��
                    call five_smooth_space_in(ub, nx, ny, S_space)
                    call five_smooth_space_in(vb, nx, ny, S_space)
                    call five_smooth_space_in(zb, nx, ny, S_space)  !  ֻƽ���˵�ǰʱ��㣬n-1ʱ��㲢δƽ������
                    
                    if (Plus_minus_smooth == 0) then  ! ��������ƽ����ֻ����ƽ��
                        continue  ! �����������ִ�У��൱��python�е�pass
                        
                    else if (Plus_minus_smooth == 1) then
                        call five_smooth_space_in(ub, nx, ny, -1 * s_space)
                        call five_smooth_space_in(vb, nx, ny, -1 * s_space)
                        call five_smooth_space_in(zb, nx, ny, -1 * s_space)  ! ��ƽ��
                    endif
                else 
                    if (mod(j * dt, 3600) == 0) then  !  ÿ1h�߽�ƽ��һ��
                        call nine_smooth_space_out(ub, nx, ny, S_space)
                        call nine_smooth_space_out(vb, nx, ny, S_space)
                        call nine_smooth_space_out(zb, nx, ny, S_space)
                        
                        if (Plus_minus_smooth == 0) then  ! ��������ƽ����ֻ����ƽ��
                            continue  ! �����������ִ�У��൱��python�е�pass
                            
                        else if (Plus_minus_smooth == 1) then  ! ������ƽ��
                            call nine_smooth_space_out(ub, nx, ny, -1 * S_space)  ! ��ƽ��
                            call nine_smooth_space_out(vb, nx, ny, -1 * S_space)
                            call nine_smooth_space_out(zb, nx, ny, -1 * S_space)
                        endif
                    endif
                endif
                
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�ж��Ƿ������ֹ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (i == int .and. j == remainder) then
                    print*, j, 'stop'
                    write(19), ((ub(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((vb(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((zb(p, q), p= 1, nx), q = 1, ny)
                    close(19)
                    stop  !  ��ֹ���� 
                endif
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���Ԥ�����!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (mod((i - 1) * 12 * 3600 + j * dt, output_interval) == 0) then
                    print*, j
                    write(19), ((ub(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((vb(p, q), p= 1, nx), q = 1, ny)
                    write(19), ((zb(p, q), p= 1, nx), q = 1, ny)
                endif
            enddo
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�ر�����ļ�!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        close(19)
    endif
END
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!����ֵ�����ӳ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!����ѭ�������Ѷ��������ʱ����������!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE transmit(ua, va, za, ub, vb, zb, nx, ny)
! ��ub, vb, zb��ֵ�ֱ𴫵ݸ�ua, va, za
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
    
    

    
    
    