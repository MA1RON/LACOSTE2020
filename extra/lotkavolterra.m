clc; clear; close all;
%% dati
aa = 1; bb = .01; cc = 1; dd = .02;

fx = @(x,y) aa*x-bb*x.*y;
fy = @(x,y) -cc*y+dd*x.*y;

x0 = 20; y0 = 20;

%% time
tend = 15;
dt = tend/1000;
time = (0:dt:tend)';
mm = length(time);

%% FE
xx = ones(mm,1)*x0; yy = ones(mm,1)*y0;
for inst = 2:mm
    xx(inst) = xx(inst-1) + dt*fx(xx(inst-1),yy(inst-1));
    yy(inst) = yy(inst-1) + dt*fy(xx(inst-1),yy(inst-1));
end

figure(1)
plot(time,xx,'b-.','linewidth',2,'DisplayName','Prede - FE')
hold on
plot(time,yy,'r-.','linewidth',2,'DisplayName','Predatori - FE')
hold on

figure(2)
plot(xx,yy,'b-.','linewidth',2,'DisplayName','FE')
hold on

%% BE
xx = ones(mm,1)*x0; yy = ones(mm,1)*y0;
optimoptions.Display = 'off';
for inst = 2:mm
    out = fsolve(@(in) [in(1)-xx(inst-1);in(2)-yy(inst-1)]-dt*[fx(in(1),in(2));fy(in(1),in(2))],...
        [xx(inst-1);yy(inst-1)],optimoptions);
    xx(inst) = out(1);
    yy(inst) = out(2);
end

figure(1)
plot(time,xx,'b--','linewidth',2,'DisplayName','Prede - BE')
hold on
plot(time,yy,'r--','linewidth',2,'DisplayName','Predatori - BE')

figure(2)
plot(xx,yy,'r-.','linewidth',2,'DisplayName','BE')

%% CN
xx = ones(mm,1)*x0; yy = ones(mm,1)*y0;
optimoptions.Display = 'off';
for inst = 2:mm
    out = fsolve(@(in) [in(1)-xx(inst-1);in(2)-yy(inst-1)]...
        - dt/2*[fx(in(1),in(2));fy(in(1),in(2))]...
        - dt/2*[fx(xx(inst-1),yy(inst-1));fy(xx(inst-1),yy(inst-1))],...
        [xx(inst-1);yy(inst-1)],optimoptions);
    xx(inst) = out(1);
    yy(inst) = out(2);
end

figure(1)
plot(time,xx,'c-.','linewidth',2,'DisplayName','Prede - CN')
hold on
plot(time,yy,'m-.','linewidth',2,'DisplayName','Predatori - CN')

figure(2)
plot(xx,yy,'m-.','linewidth',2,'DisplayName','CN')

%% ---
figure(1)
legend('location','northeastoutside')
set(gca,'fontsize',18)
grid on
box on
xlabel('Tempo')
ylabel('Numero di individui')

figure(2)
legend('location','northeastoutside')
set(gca,'fontsize',18)
grid on
box on
xlabel('# Prede')
ylabel('# Predatori')