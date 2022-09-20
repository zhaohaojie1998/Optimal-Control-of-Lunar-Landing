close all

% ######### 绘图 ######### %

%--------------------------------
 reset(groot); % 重置背景色
figure(1)
plot(data.t,data.tgo,'k--',LineWidth=2)
legend('待飞时间')

axis square
grid on

title('待飞时间估计')
xlabel('Time (s)') 
ylabel('Time to go (s)')
%--------------------------------


%--------------------------------
figure(2)
plot(data.crossrange/1e3,data.downrange/1e3,'k--',LineWidth=2)
hold on
plot(data.crossrange(1)/1e3,data.downrange(1)/1e3,'ko',LineWidth=1,MarkerSize=10)
hold on
plot(0,0,'p',color='k',LineWidth=1,MarkerSize=10)
legend('相对偏差','初始估计','目标落点')

axis square
grid on

title('纵向-横向偏差估计')
xlabel('Cross-range (km)') 
ylabel('Down-range (km)')
%--------------------------------


%--------------------------------
figure(3)
plot(data.t,data.altitude/1e3,'k-.',LineWidth=2)
legend('高度')

grid on

title('高度-时间曲线')
xlabel('Time (s)') 
ylabel('Altitude (km)')
%--------------------------------


%--------------------------------
figure(4)
plot(rad2deg(data.longitude),rad2deg(data.latitude),'k-.',LineWidth=2)
hold on
plot(rad2deg(f.Lon),rad2deg(f.Lat),'p',color='k',LineWidth=1,MarkerSize=10)
legend('星下点轨迹','目标着陆点')

axis square
%axis equal
grid on

title('星下点轨迹')
xlabel('Longitude (deg)') 
ylabel('Latitude (deg)')
%--------------------------------


%--------------------------------
figure(12)
plot3(data.xyz(:,1)/1e3,data.xyz(:,2)/1e3,data.xyz(:,3)/1e3,'k-.',LineWidth=2)
hold on 
plot3(0,0,0,'p',color='k',LineWidth=1,MarkerSize=10)
legend('飞 行 轨 迹','目标着陆点')

set(gca,'XDir','reverse')%对x方向反转
set(gca,'ZDir','reverse')%对z方向反转
% axis equal
axis square
grid on

title('目标北东地坐标系着陆飞行轨迹')
xlabel('North (km)') 
ylabel('East (km)')
zlabel('Ground (km)')
%--------------------------------



%--------------------------------
figure(11)
plot(data.t,data.vx/1e3,'-',data.t,data.vy/1e3,'--',data.t,data.vz/1e3,'-.',LineWidth=1)
legend('V_x','V_y','V_z')

grid on

title('月心惯性系速度时间曲线')
xlabel('Time (s)') 
ylabel('Velocity (km/s)')
%--------------------------------


%--------------------------------
figure(10)
plot(data.t,data.m,'r',LineWidth=1)
legend('m')

grid on

title('质量曲线')
xlabel('Time (s)') 
ylabel('Mass (kg)')
%--------------------------------


%--------------------------------
figure(5)
plot(data.t,data.T)
legend('T')

grid on

title('推力大小')
xlabel('Time (s)') 
ylabel('Thrust (N)')
%--------------------------------


%--------------------------------
figure(6)
plot(data.t,rad2deg(data.iT(:,1)),'-',data.t,rad2deg(data.iT(:,2)),'--',LineWidth=1)
legend('T 与 Ox 轴 夹 角','T 与 xy 平面夹角')

grid on

title('惯性系推力方向')
xlabel('Time (s)') 
ylabel('Thrust Direction (deg)')
%--------------------------------


%--------------------------------
figure(7)
plot(data.t,data.k*100,'g-.')
legend('k')

grid on

title('油门大小')
xlabel('Time (s)') 
ylabel('Throttle (%)')
%--------------------------------


%--------------------------------
figure(8)
plot(data.t,data.u(:,1),'-',data.t,data.u(:,2),'--',data.t,data.u(:,3),'-.')
legend('u_x','u_y','u_z')

grid on

title('控制量（实际推力加速度）')
xlabel('Time (s)') 
ylabel('Control (m/s^2)')
%--------------------------------


%--------------------------------
figure(9)
plot(data.t,data.aT(:,1),'-',data.t,data.aT(:,2),'--',data.t,data.aT(:,3),'-.',LineWidth=1)
legend('max a_x','max a_y','max a_z')

grid on

title('最大推力加速度（推力方向）')
xlabel('Time (s)') 
ylabel('Max Thrust Acceleration (m/s^2)')
%--------------------------------



%--------------------------------
figure(13)
colordef none;  %2D/3D图窗背景透明

plot3(o.orbit(:,1)/1e3,o.orbit(:,2)/1e3,o.orbit(:,3)/1e3,color ='#ADFF2F',lines='-.',LineWidth=2)
hold on 
plot3(data.rx/1e3,data.ry/1e3,data.rz/1e3,color = '#FFD700',LineWidth=2)
hold on 
plot3(f.rx/1e3,f.ry/1e3,f.rz/1e3,'p',color='r',LineWidth=2,MarkerSize=5)
hold on

% 绘制月球
PIC = imread('LunarSurface.jpg');
[x,y,z] = sphere(50);                      % 创建球体
x = Radius*x*1e-3;
y = Radius*y*1e-3;
z = Radius*z*1e-3;
s = surface(x,y,z);                          % 绘制球面
s.FaceColor = 'texturemap';         % 使用纹理贴图
s.CData = PIC;                % 将颜色数据设置为地形数据
s.EdgeColor = 'none';                   % 删除 edges
s.FaceLighting = 'gouraud';         % 曲面的首选照明百分比
s.SpecularStrength = 0.5;             % 更改反射光的强度
light('Position',[2 0 0])                % 添加灯光
view([60,10])                                % 设置视角
hold on

% 月球框架
% plot3(0,0,0,'o',color='#778899',LineWidth=1,MarkerSize=8)
% hold on
% plt.theta = -pi:pi/50:pi;
% plt.x = Radius*cos(plt.theta);
% plt.z = Radius*sin(plt.theta);
% plt.y = 0.*plt.theta;
% for i = 15:15:360
% [plt.x,plt.y,plt.z] = plotMoon(plt.x,plt.y,plt.z, i);
% plot3(plt.x/1e3,plt.y/1e3,plt.z/1e3,color='#778899',LineWidth=1,LineStyle=':')
% hold on
% end
% for i = 0:15:90
% plt.z = Radius*sin(deg2rad(i))*ones(length(plt.theta));
% plt.x = Radius*cos(deg2rad(i))*cos(plt.theta);
% plt.y = Radius*cos(deg2rad(i))*sin(plt.theta);
% plot3(plt.x/1e3,plt.y/1e3,plt.z/1e3,color='#778899',LineWidth=1,LineStyle=':')
% hold on
% plot3(plt.x/1e3,plt.y/1e3,-plt.z/1e3,color='#778899',LineWidth=1,LineStyle=':')
% hold on
% end

set(gcf,'color','#000033'); %设置显示界面背景
legend('初 始 轨 道','着 陆 轨 迹','目标着陆点')
axis square
% axis square
grid on
title('月球软着陆轨道')
xlabel('X (km)') 
ylabel('Y (km)')
zlabel('Z (km)')
%--------------------------------







function [x,y,z] = plotMoon(x,y,z,lambda)
lambda = deg2rad(lambda);
    A = [ cos(lambda), sin(lambda), 0; 
              -sin(lambda), cos(lambda), 0;
              0,0,1;];
    xyz = A*[x;y;z];
    x = xyz(1,:);
    y = xyz(2,:);
    z = xyz(3,:);

end