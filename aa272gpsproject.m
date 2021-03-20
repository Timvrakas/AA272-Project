
clear
close all

format long

%%

%import


%s = readtable("ephem.csv");
%t = readtable("gnss_log.csv");

% s = read_rinex_nav('idtd0710.21n');

%Any CORS nav files
s = read_rinex_nav('guam2010.00n');
%RINEX 2.11 (2 variants, use 0 or 1)
% t = read_rinex_obs('idtd0710.21o',0);
%RINEX 2.10
t = read_rinex_obs_legacy('guam2010.00o');
%RINEX 3.03
%t = read_phone_rinex_obs('Pixel4_GnssLog.20o');

%%

% [Weekmin tmin] = UT1toGPStime(06,08,20,22,52,28);
% [Weekmax tmax] = UT1toGPStime(06,08,20,23,21,32.4334262);

for n = 1:length(t.tr)
    
    Nw = t.GPSWeek(n);
    
    [x_ecef(n), y_ecef(n), z_ecef(n), b(n)] = calculateSatPos(s,Nw,t.prn(n),t.tr(n));
    
end

t.X = x_ecef';
t.Y = y_ecef';
t.Z = z_ecef';
t.B = b';

figure
plot3(t.X,t.Y,t.Z,'.')
xlabel("X (km)")
ylabel("Y (km)")
zlabel("Z (km)")
title("Satellite positions")
axis  equal;
view (3);


%%

n = 1;
p = 1;
e = 1*(10^-8);
change_lim = 10000000;
y(:,p) = [0;0;0;0];
leng = length(t.X);
r = 0;
IF = table();

% nmin = find(t.tr(t.GPSWeek == Weekmin) >= tmin,1,'first');
% nmax = find(t.tr(t.GPSWeek == Weekmax) <= tmax,1,'last');

time(p) = t.tr(n);

%%while(n<nmax)
while(n<leng)
    
    [X,Y,Z,B,P,n,time(p+1),svid,IFnew] = findNextTimeIF(n,t);
    
    see = length(X);
    
    if(length(X) >= 4)
        y_temp = NewtonRhapson(y(:,p),X,Y,Z,B,P,e);
        if(norm(y(:,p) - y_temp)< change_lim)||(p == 1)
            y(:,p+1) = y_temp;
            if(IFnew.Svid(1)~=0)
                IFnew.x = y(1,p+1)*ones(size(IFnew.X));
                IFnew.y = y(2,p+1)*ones(size(IFnew.X));
                IFnew.z = y(3,p+1)*ones(size(IFnew.X));
                IF = [IF;IFnew];
            end
            p = p+1;

        else
            r = r+1;
            reject(1:4,r) = y_temp;
            reject(5,r) = see;
            reject(6,r) = time(p+1);
        end
    end
    
end

%%

IF = mapIono(IF);

%%

figure()    

scatter(rad2deg(IF.long_bot),rad2deg(IF.lat_bot),40,IF.A,'filled')
hold on
colorbar

colormap jet;
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = "A (m*HZ^2)";

% Load and plot MATLAB built-in Earth topography data
load('topo.mat', 'topo');
topoplot = [topo(:, 181:360), topo(:, 1:180)];
contour(-180:179, -90:89, topoplot, [0, 0], 'black');

axis equal
grid on
xlim([-180, 180]);
ylim([-90, 90]);
xlabel('Longitude [\circ]');
ylabel('Latitude [\circ]');

figure()

for n = 1:length(IF.X)

    plot([rad2deg(IF.long_bot(n)) rad2deg(IF.long_top(n))],[rad2deg(IF.lat_bot(n)) rad2deg(IF.lat_top(n))],'r')
    hold on
    
end

% Load and plot MATLAB built-in Earth topography data
load('topo.mat', 'topo');
topoplot = [topo(:, 181:360), topo(:, 1:180)];
contour(-180:179, -90:89, topoplot, [0, 0], 'black');

axis equal
grid on
xlim([-180, 180]);
ylim([-90, 90]);
xlabel('Longitude [\circ]');
ylabel('Latitude [\circ]');

Re = 6378137; %m
[xE , yE, zE] = ellipsoid(0, 0, 0, Re , Re, Re, 20);

figure(10)
surface(xE, yE , zE ,'FaceColor','blue','EdgeColor','black');
hold on

for n = 1:length(IF.A)

    plot3([IF.x_bot_ECEF(n) IF.x_top_ECEF(n)],[IF.y_bot_ECEF(n) IF.y_top_ECEF(n)],[IF.z_bot_ECEF(n) IF.z_top_ECEF(n)],'r')
    hold on
end


axis  equal;
view (3);
grid on;
xlabel("X (km)")
ylabel("Y (km)")
zlabel("Z (km)")


time_rel = (time(2:end) - time(2))/(10^3);

figure()
scatter3(y(1,2:end)/1000,y(2,2:end)/1000,y(3,2:end)/1000,40,time_rel,'filled')
colorbar
xlabel("X (km)")
ylabel("Y (km)")
zlabel("Z (km)")
title("New Trajectory")
axis  equal;
view (3);

colormap jet;
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = "Time (s)";

writetable(IF,"IF_output.csv")

%IF

function IF = mapIono(IF)

    Re = 6378137; %m
    R1 = 50*1000 + Re; %m
    R2 = 1000*1000 + Re; %m

    for n = 1:length(IF.X)
        [el(n,1),az(n,1)] = getAngles(IF.X(n),IF.Y(n),IF.Z(n),[IF.x(n) IF.y(n) IF.z(n)]); 
    end

    a = pi/2 - el - asin(Re*sin(el + pi/2)/R1);
    b = pi/2 - el - asin(Re*sin(el + pi/2)/R2);
    d1 = R1*sin(a)./sin(el + pi/2);
    d2 = R2*sin(b)./sin(el + pi/2);
%    A_norm = IF.A./(d2-d1);
%     lat0 = asin(IF.z./sqrt((IF.x).^2 + (IF.y).^2 + (IF.z).^2));
%     long0 = atan2(IF.y,IF.x);
    d = sqrt((IF.X-IF.x).^2 + (IF.Y-IF.y).^2 + (IF.Z-IF.z).^2);
    
%     ra = Re*a;
%     dxa = ra.*cos(az);
%     dya = ra.*sin(az);
%     dlata = 2*pi*dya/Re;
%     dlonga = 2*pi*dxa./(Re*cos(lat0));
%     lata = lat0 + dlata;
%     longa = long0 + dlonga;
    x_bot_ECEF = (d1./d).*(IF.X-IF.x) + IF.x;
    y_bot_ECEF = (d1./d).*(IF.Y-IF.y) + IF.y;
    z_bot_ECEF = (d1./d).*(IF.Z-IF.z) + IF.z;
    lata = asin(z_bot_ECEF./sqrt(x_bot_ECEF.^2 + y_bot_ECEF.^2 + z_bot_ECEF.^2));
    longa = atan2(y_bot_ECEF,x_bot_ECEF);
    
%     rb = Re*b;
%     dxb = rb.*cos(az);
%     dyb = rb.*sin(az);
%     dlatb = 2*pi*dyb/Re;
%     dlongb = 2*pi*dxb./(Re*cos(lat0));
%     latb = lat0 + dlatb;
%     longb = long0 + dlongb;
    x_top_ECEF = (d2./d).*(IF.X-IF.x) + IF.x;
    y_top_ECEF = (d2./d).*(IF.Y-IF.y) + IF.y;
    z_top_ECEF = (d2./d).*(IF.Z-IF.z) + IF.z;
    latb = asin(z_top_ECEF./sqrt(x_top_ECEF.^2 + y_top_ECEF.^2 + z_top_ECEF.^2));
    longb = atan2(y_top_ECEF,x_top_ECEF);
    
%    IF.A_norm = A_norm;
    IF.lat_bot = lata;
    IF.long_bot = longa;
    IF.lat_top = latb;
    IF.long_top = longb;
    IF.x_bot_ECEF = x_bot_ECEF;
    IF.y_bot_ECEF = y_bot_ECEF;
    IF.z_bot_ECEF = z_bot_ECEF;
    IF.x_top_ECEF = x_top_ECEF;
    IF.y_top_ECEF = y_top_ECEF;
    IF.z_top_ECEF = z_top_ECEF;
    
end

%Orbit functions

function [el, az] = getAngles(X,Y,Z,y)
    r_ecef = [X-y(1), Y- y(2), Z-y(3)]';
    lat = asin(y(3)/norm(y(1:3)));
    long = atan2(y(2),y(1));
    
    E = [-sin(long) cos(long) 0];
    N = [-sin(lat)*cos(long) -sin(lat)*sin(long) cos(lat)];
    U = [cos(lat)*cos(long) cos(lat)*sin(long) sin(lat)];
    R = [E; N; U];
    
    r_enu = R*r_ecef;
    
    az = atan2(r_enu(1),r_enu(2));
    el = atan2(r_enu(3),norm(r_enu(1:2)));
    
end

function E_end = NewtonRaphsonE(M,epsilon,e)
    
    E = zeros(1,51);

    if(e <.5)
        E(1) = pi;
    else
        E(1) = M;
    end
    
    for i = 1:50
        d = (E(i)-e*sin(E(i))-M)/(1-e*cos(E(i)));
        E(i+1) = E(i) - d;
        
        if abs(d) < epsilon
            E(i+1) = E(i);
        end
    end
    
    E_end = E(end);

end

function [x_ecef y_ecef z_ecef b] = calculateSatPos(s,Nw,prn,t_in)

    mu = 3.986005*(10^14); %m^3/s^2
    Om_dot_e = 7.2921151467*(10^-5); %rad/s
    c = 299792458; %m/s
    F = -4.442807633*(10^-10);
    
    
    tdiff = (Nw - s.GPSWeekToc)*604800 + (t_in - s.Toc);
    rows_available = find(tdiff >= 0);
    rows_sat = find(s.PRN(rows_available) == prn);
    if length(rows_sat) < 1
        message = "Error"
        rows_sat = find(s.PRN == prn);
        id = rows_sat(1);
    else
        id = rows_sat(end);   
    end

    a = (s.sqrtA(id))^2;
    n = sqrt(mu/(a^3));
    e = s.Eccentricity(id);

    tk = (Nw - s.GPSWeek(id))*604800 + (t_in - s.Toe(id));
    Mk = s.M0(id) + n.*tk;
    Ek = NewtonRaphsonE(Mk,10^-9,e);
    s_vk = (sqrt(1-e^2)*sin(Ek))/(1-e*cos(Ek));
    c_vk = (cos(Ek) - e)/(1-e*cos(Ek));
    vk = atan2(s_vk,c_vk);
    phik = vk + s.omega(id);

    uk = phik;
    for i = 1:50

        dphik = s.Cus(id)*sin(2*uk) + s.Cuc(id)*cos(2*uk);
        uk = phik + dphik;

        store(1:2,i) = [dphik;uk];

    end

    drk = s.Crs(id)*sin(2*phik) + s.Crc(id)*cos(2*phik);
    dik = s.Cis(id)*sin(2*phik) + s.Cic(id)*cos(2*phik);
    Omk = s.Omega0(id) - Om_dot_e*t_in + s.OmegaDot(id)*tk;
    rk = a*(1 - e*cos(Ek)) + drk;
    ik = s.Io(id) + s.IDOT(id)*tk + dik;
    xp = rk*cos(uk);
    yp = rk*sin(uk);
    x_ecef = xp*cos(Omk) - yp*cos(ik)*sin(Omk);
    y_ecef = xp*sin(Omk) + yp*cos(ik)*cos(Omk);
    z_ecef = yp*sin(ik);
    b = c*(s.af0(id)+s.af1(id)*(tdiff(id))+s.af2(id)*((tdiff(id))^2) + F*e*sqrt(a)*sin(Ek));

end


%Positioning functions

function [X,Y,Z,B,P,n_f,time,ID,IF] = findNextTimeIF(n_i,t)

    f1_2 = (1575.42*(10^6))^2;
    f2_2 = (1227.60*(10^6))^2;
    f5_2 = (1176.45*(10^6))^2;
    a1 = f1_2/(f1_2 - f2_2);
    b1 = f2_2/(f1_2 - f2_2);
    c1 = f1_2*f2_2/(f2_2 - f1_2);
    a2 = f1_2/(f1_2 - f5_2);
    b2 = f5_2/(f1_2 - f5_2);
    c2 = f1_2*f5_2/(f5_2 - f1_2);
    a3 = f2_2/(f2_2 - f5_2);
    b3 = f5_2/(f2_2 - f5_2);
    c3 = f2_2*f5_2/(f5_2 - f2_2);
    
    IF = table();
    time = t.tr(n_i);
    leng = length(t.X);
    n = n_i;
    t_check = time;
    s = 1;
    m = 1;
    ds = 100;
    
    while (t_check == time)&&(n<leng)
        
        if(~isnan(t.P_L1(n)) && ~isnan(t.P_L2(n)))
            X(s) = t.X(n);
            Y(s) = t.Y(n);
            Z(s) = t.Z(n);
            B(s) = t.B(n);
            P(s) = a1*t.P_L1(n) - b1*t.P_L2(n);
            svid(s) = t.prn(n);
            AIF(m) = c1*(t.P_L1(n) - t.P_L2(n));
            PseuIF(m) = P(s);
            timeIF(m) = t.tr(n);
            svidIF(m)= svid(s); 
            XIF(m) = X(s);
            YIF(m) = Y(s);
            ZIF(m) = Z(s);
            m = m+1;
        elseif(~isnan(t.P_L1(n)) && ~isnan(t.P_L5(n)))
            X(s) = t.X(n);
            Y(s) = t.Y(n);
            Z(s) = t.Z(n);
            B(s) = t.B(n);
            P(s) = a2*t.P_L1(n) - b2*t.P_L5(n);
            svid(s) = t.prn(n);
            AIF(m) = c2*(t.P_L1(n) - t.P_L5(n));
            PseuIF(m) = P(s);
            timeIF(m) = t.tr(n);
            svidIF(m)= svid(s); 
            XIF(m) = X(s);
            YIF(m) = Y(s);
            ZIF(m) = Z(s);
            m = m+1;
        elseif(~isnan(t.P_L2(n)) && ~isnan(t.P_L5(n)))
            X(s) = t.X(n);
            Y(s) = t.Y(n);
            Z(s) = t.Z(n);
            B(s) = t.B(n);
            P(s) = a3*t.P_L2(n) - b3*t.P_L5(n);
            svid(s) = t.prn(n);
            AIF(m) = c3*(t.P_L2(n) - t.P_L5(n));
            PseuIF(m) = P(s);
            timeIF(m) = t.tr(n);
            svidIF(m)= svid(s); 
            XIF(m) = X(s);
            YIF(m) = Y(s);
            ZIF(m) = Z(s);
            m = m+1;
        else
            if ~isnan(t.P_L1(n))
                X(s) = t.X(n);
                Y(s) = t.Y(n);
                Z(s) = t.Z(n);
                B(s) = t.B(n);
                P(s) = t.P_L1(n);
                svid(s) = t.prn(n);
            elseif ~isnan(t.P_L2(n))
                X(s) = t.X(n);
                Y(s) = t.Y(n);
                Z(s) = t.Z(n);
                B(s) = t.B(n);
                P(s) = t.P_L2(n);
                svid(s) = t.prn(n);
            elseif ~isnan(t.P_L5(n))
                X(s) = t.X(n);
                Y(s) = t.Y(n);
                Z(s) = t.Z(n);
                B(s) = t.B(n);
                P(s) = t.P_L5(n);
                svid(s) = t.prn(n);
            else
                s = s-1
            end  
        end
        
        s = s+1;
        
        n = n + 1;
        t_check = t.tr(n);
    end
    
    ID = svid;
    n_f = n;
    
    if(m>1)
        IF.Svid = svidIF';
        IF.time = timeIF';
        IF.A = AIF';
        IF.P = PseuIF';
        IF.X = XIF';
        IF.Y = YIF';
        IF.Z = ZIF';
    else
        IF.Svid = 0;
        IF.time = 0;
        IF.A = 0;
        IF.P = 0;
        IF.X = 0;
        IF.Y = 0;
        IF.Z = 0;
    end
    
      
end

function y_end = NewtonRhapson(y0,X,Y,Z,B,P,e)

    y(:,1) = y0;
    
    dy_norm = 1000;
    n = 1;
    
    while(abs(dy_norm) > e)&&(n<100)
        
        dy = positionUpdate(y(:,n),X,Y,Z,B,P);
        
        y(:,n+1) = y(:,n) + dy;
        dy_norm = norm(dy);
        n = n + 1;
        
    end
    
    y_end = y(:,n);
    
end

function dy = positionUpdate(y0,X,Y,Z,B,p)
    
    x = y0(1);
    y = y0(2);
    z = y0(3);
    bu = y0(4);
    
    p_ex = expectedPseudorange(x,y,z,bu,X,Y,Z,B);
    G = constructGPos(x,y,z,X,Y,Z);
    
    dp = p' - p_ex;
    
    dy = inv(G'*G)*G'*dp;
    
end

function p = expectedPseudorange(x,y,z,bu,X,Y,Z,B)

    p = zeros(length(X),1);
    
    for i = 1:length(X)
        
        xdiff = [X(i)-x;
            Y(i) - y;
            Z(i) - z];
        
        p(i) = norm(xdiff) + bu - B(i);
        
    end

end

function G = constructGPos(x,y,z,X,Y,Z)

    G = zeros(length(X),4);
    
    for i = 1:length(X)
        
        xdiff = [X(i)-x;
            Y(i) - y;
            Z(i) - z];
        p = norm(xdiff);
        l = xdiff/p;
        
        G(i,:) = [-l' 1];
        
    end

end

%RINEX
%                  col  1:    PRN    ....... satellite PRN          
%                  col  2:    M0     ....... mean anomaly at reference time
%                  col  3:    delta_n  ..... mean motion difference
%                  col  4:    e      ....... eccentricity
%                  col  5:    sqrt(A)  ..... where A is semimajor axis
%                  col  6:    OMEGA  ....... LoAN at weekly epoch
%                  col  7:    i0     ....... inclination at reference time
%                  col  8:    omega  ....... argument of perigee
%                  col  9:    OMEGA_dot  ... rate of right ascension 
%                  col 10:    i_dot  ....... rate of inclination angle
%                  col 11:    Cuc    ....... cosine term, arg. of latitude
%                  col 12:    Cus    ....... sine term, arg. of latitude
%                  col 13:    Crc    ....... cosine term, radius
%                  col 14:    Crs    ....... sine term, radius
%                  col 15:    Cic    ....... cosine term, inclination
%                  col 16:    Cis    ....... sine term, inclination
%                  col 17:    toe    ....... time of ephemeris
%                  col 18:    IODE   ....... Issue of Data Ephemeris
%                  col 19:    GPS_wk ....... GPS week
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = read_rinex_nav( filename )
    s = table();
    fid = fopen(filename);
    if fid == -1
        errordlg(['The file ''' filename ''' does not exist.']);
        return;
    end
    % skip through header
    end_of_header = 0;
    while end_of_header == 0
        current_line = fgetl(fid);
        if contains(current_line,'END OF HEADER')
            end_of_header=1;
        end
    end
    j = 0;
    while feof(fid) ~= 1
        j = j+1;

        current_line = fgetl(fid);
        % parse epoch line (ignores SV clock bias, drift, and drift rate)
        PRN(j) = str2double(current_line(1:2));
        Y(j) = str2double(current_line(3:5));
        M(j) = str2double(current_line(6:8));
        D(j) = str2double(current_line(9:11));
        H(j) = str2double(current_line(12:14));
        min(j) = str2double(current_line(15:17));
        sec(j) = str2double(current_line(18:22));
        af0(j) = str2double(current_line(23:37))*(10^str2double(current_line(39:41)));
        af1(j) = str2double(current_line(42:56))*(10^str2double(current_line(58:60)));
        af2(j) = str2double(current_line(61:75))*(10^str2double(current_line(77:79)));
        
        % Broadcast orbit line 1
        current_line = fgetl(fid);
        IODE(j) = str2double(current_line(4:18))*(10^str2double(current_line(20:22)));
        Crs(j) = str2double(current_line(23:37))*(10^str2double(current_line(39:41)));
        deltaN(j) = str2double(current_line(42:56))*(10^str2double(current_line(58:60)));
        M0(j) = str2double(current_line(61:75))*(10^str2double(current_line(77:79)));
        % Broadcast orbit line 2
        current_line = fgetl(fid);
        Cuc(j) = str2double(current_line(4:18))*(10^str2double(current_line(20:22)));
        Eccentricity(j) = str2double(current_line(23:37))*(10^str2double(current_line(39:41)));
        Cus(j) = str2double(current_line(42:56))*(10^str2double(current_line(58:60)));
        sqrtA(j) = str2double(current_line(61:75))*(10^str2double(current_line(77:79)));
        % Broadcast orbit line 3
        current_line = fgetl(fid);
        Toe(j) = str2double(current_line(4:18))*(10^str2double(current_line(20:22)));
        Cic(j) = str2double(current_line(23:37))*(10^str2double(current_line(39:41)));
        Omega0(j) = str2double(current_line(42:56))*(10^str2double(current_line(58:60)));
        Cis(j) = str2double(current_line(61:75))*(10^str2double(current_line(77:79)));
        % Broadcast orbit line 4
        current_line = fgetl(fid);
        Io(j) = str2double(current_line(4:18))*(10^str2double(current_line(20:22)));
        Crc(j) = str2double(current_line(23:37))*(10^str2double(current_line(39:41)));
        omega(j) = str2double(current_line(42:56))*(10^str2double(current_line(58:60)));
        OmegaDot(j) = str2double(current_line(61:75))*(10^str2double(current_line(77:79)));
        % Broadcast orbit line 5
        current_line = fgetl(fid);
        IDOT(j) = str2double(current_line(4:18))*(10^str2double(current_line(20:22)));
        L2_codes(j) = str2double(current_line(23:37))*(10^str2double(current_line(39:41)));
        GPSWeek(j) = str2double(current_line(42:56))*(10^str2double(current_line(58:60)));
        L2_dataflag(j) = str2double(current_line(61:75))*(10^str2double(current_line(77:79)));
        % Broadcast orbit line 6
        current_line = fgetl(fid);
        SV_acc(j) = str2double(current_line(4:18))*(10^str2double(current_line(20:22)));
        SV_health(j) = str2double(current_line(23:37))*(10^str2double(current_line(39:41)));
        TGD(j) = str2double(current_line(42:56))*(10^str2double(current_line(58:60)));
        IODC(j) = str2double(current_line(61:75))*(10^str2double(current_line(77:79)));
        % Broadcast orbit line 7
        current_line = fgetl(fid);
        msg_trans_t(j) = str2double(current_line(4:18))*(10^str2double(current_line(20:22)));
        if (length(current_line) >= 41)
            fit_int(j) = str2double(current_line(23:37))*(10^str2double(current_line(39:41)));
        else
            fit_int(j) = 0;
        end
    end
    
    s.PRN = PRN';
    s.year = Y';
    s.M = M';
    s.D = D';
    s.H = H';
    s.min = min';
    s.sec = sec';
    [GPSWeekToc Toc] = UT1toGPStime(D,M,Y,H,min,sec);
    s.GPSWeekToc = GPSWeekToc;
    s.Toc = Toc;
    s.af0 = af0';
    s.af1 = af1';
    s.af2 = af2';
    s.IODE = IODE';
    s.Crs = Crs';
    s.deltaN = deltaN';
    s.M0 = M0';
    s.Cuc = Cus';
    s.Eccentricity = Eccentricity';
    s.Cus = Cus';
    s.sqrtA = sqrtA';
    s.Toe = Toe';
    s.Cic = Cic';
    s.Omega0 = Omega0';
    s.Cis = Cis';
    s.Io = Io';
    s.Crc = Crc';
    s.omega = omega';
    s.OmegaDot = OmegaDot';
    s.IDOT = IDOT';
    s.L2_codes = L2_codes';
    s.GPSWeek = GPSWeek';
    s.L2_dataflag = L2_dataflag';
    s.SV_acc = SV_acc';
    s.SV_health = SV_health';
    s.TGD = TGD';
    s.IODC = IODC';
    s.msg_trans_t = msg_trans_t';
    s.fit_int = fit_int';

end
function t = read_rinex_obs( filename,extra )
    t = table();
    fid = fopen(filename);
    if fid == -1
        errordlg(['The file ''' filename ''' does not exist.']);
        return;
    end
    % skip through header
    end_of_header = 0;
    while end_of_header == 0
        current_line = fgetl(fid);
        if contains(current_line,'END OF HEADER')
            end_of_header=1;
        end
    end
    j = 0;
    n = 1;
    while feof(fid) ~= 1
        j = j+1;

        current_line = fgetl(fid);
        
        if (length(current_line) < 33)
            current_line = fgetl(fid);
            current_line = fgetl(fid);
            current_line = fgetl(fid);
            current_line = fgetl(fid);
            current_line = fgetl(fid);
        end
        
        if (length(current_line) < 33)
            F = '?'
        end
        
        Y_t = str2double(current_line(1:3));
        M_t = str2double(current_line(4:6));
        D_t = str2double(current_line(7:9));
        H_t = str2double(current_line(10:12));
        min_t = str2double(current_line(13:15));
        sec_t = str2double(current_line(16:26));
        epoch_flag_t = str2double(current_line(27:29));
        num_sat = str2double(current_line(30:32));
        if(num_sat <=12)
            prns = current_line(33:(32+3*num_sat));
            prn_list = strsplit(prns, 'G');
            prn_array = str2double(prn_list(2:(num_sat+1)));
        else
            prns = current_line(33:(32+3*12));
            prn_list = strsplit(prns, 'G');
            prn_array = str2double(prn_list(2:(12+1)));
            current_line = fgetl(fid);
            prns = current_line(33:(32+3*(num_sat-12)));
            prn_list = strsplit(prns, 'G');
            prn_array = [prn_array str2double(prn_list(2:(num_sat-11)))];  
        end
        
        
        for i = 1:num_sat
            % Broadcast orbit line 1
            current_line = fgetl(fid);
%             if (length(current_line) >= 14) L1(n) = str2double(current_line(1:14)); else L1(n) = NaN; end
%             if (length(current_line) >= 15) L1_LLI(n) = str2double(current_line(15)); else L1_LLI(n) = NaN; end
%             if (length(current_line) >= 16) L1_str(n) = str2double(current_line(16)); else L1_str(n) = NaN; end
%             if (length(current_line) >= 30) L2(n) = str2double(current_line(17:30)); else L2(n) = NaN; end
%             if (length(current_line) >= 31) L2_LLI(n) = str2double(current_line(31)); else L2_LLI(n) = NaN; end
%             if (length(current_line) >= 32) L2_str(n) = str2double(current_line(32)); else L2_str(n) = NaN; end
%             if (length(current_line) >= 46) L5(n) = str2double(current_line(33:46)); else L5(n) = NaN; end
%             if (length(current_line) >= 47) L5_LLI(n) = str2double(current_line(47)); else L5_LLI(n) = NaN; end
%             if (length(current_line) >= 48) L5_str(n) = str2double(current_line(48)); else L5_str(n) = NaN; end
            if (length(current_line) >= 62) C1(n) = str2double(current_line(49:62)); else C1(n) = NaN; end
            if (length(current_line) >= 63) C1_LLI(n) = str2double(current_line(63)); else C1_LLI(n) = NaN; end
            if (length(current_line) >= 64) C1_str(n) = str2double(current_line(64)); else C1_str(n) = NaN; end
            if (length(current_line) >= 78) P1(n) = str2double(current_line(65:78)); else P1(n) = NaN; end
            if (length(current_line) >= 79) P1_LLI(n) = str2double(current_line(79)); else P1_LLI(n) = NaN; end
            if (length(current_line) >= 80) P1_str(n) = str2double(current_line(80)); else P1_str(n) = NaN; end
            % Broadcast orbit line 2
            current_line = fgetl(fid);
            if (length(current_line) >= 14) C2(n) = str2double(current_line(1:14)); else C2(n) = NaN; end
            if (length(current_line) >= 15) C2_LLI(n) = str2double(current_line(15)); else C2_LLI(n) = NaN; end
            if (length(current_line) >= 16) C2_str(n) = str2double(current_line(16)); else C2_str(n) = NaN; end
            if (length(current_line) >= 30) P2(n) = str2double(current_line(17:30)); else P2(n) = NaN; end
            if (length(current_line) >= 31) P2_LLI(n) = str2double(current_line(31)); else P2_LLI(n) = NaN; end
            if (length(current_line) >= 32) P2_str(n) = str2double(current_line(32)); else P2_str(n) = NaN; end
            if (length(current_line) >= 46) C5(n) = str2double(current_line(33:46)); else C5(n) = NaN; end
%             if (length(current_line) >= 47) C5_LLI(n) = str2double(current_line(47)); else C5_LLI(n) = NaN; end
%             if (length(current_line) >= 48) C5_str(n) = str2double(current_line(48)); else C5_str(n) = NaN; end
%             if (length(current_line) >= 62) S1(n) = str2double(current_line(49:62)); else S1(n) = NaN; end
%             if (length(current_line) >= 63) S1_LLI(n) = str2double(current_line(63)); else S1_LLI(n) = NaN; end
%             if (length(current_line) >= 64) S1_str(n) = str2double(current_line(64)); else S1_str(n) = NaN; end
%             if (length(current_line) >= 78) S2(n) = str2double(current_line(65:78)); else S2(n) = NaN; end
%             if (length(current_line) >= 79) S2_LLI(n) = str2double(current_line(79)); else S2_LLI(n) = NaN; end
%             if (length(current_line) >= 80) S2_str(n) = str2double(current_line(80)); else S2_str(n) = NaN; end
%             
            % Broadcast orbit line 3
            current_line = fgetl(fid);
%             if (length(current_line) >= 14) S5(n) = str2double(current_line(1:14)); else S5(n) = NaN; end
%             if (length(current_line) >= 15) S5_LLI(n) = str2double(current_line(15)); else S5_LLI(n) = NaN; end
%             if (length(current_line) >= 16) S5_str(n) = str2double(current_line(16)); else S5_str(n) = NaN; end
%             
            if(extra == 1)
                current_line = fgetl(fid);
            end
            
            Y(n) = Y_t;
            M(n) = M_t;
            D(n) = D_t;
            H(n) = H_t;
            min(n) = min_t;
            sec(n) = sec_t;
            epoch_flag(n) = epoch_flag_t;
            prn(n) = prn_array(i);
            if (~isnan(C1(n)))
                Pseu_L1(n) = C1(n);
            else
                Pseu_L1(n) = P1(n);
            end
            if (~isnan(C2(n)))
                Pseu_L2(n) = C2(n);
            else
                Pseu_L2(n) = P2(n);
            end
            Pseu_L5(n) = C5(n);
            
            n = n + 1;
            
        end
    end
    
    t.year = Y';
    t.M = M';
    t.D = D';
    t.H = H';
    t.min = min';
    t.sec = sec';
    [GPSWeek tr] = UT1toGPStime(D,M,Y,H,min,sec);
    t.GPSWeek = GPSWeek;
    t.tr = tr;
    t.ut1time = UT1tosecs(D,M,Y,H,min,sec);
    t.epoch_flag = epoch_flag';
    t.prn = prn';
    t.P_L1 = Pseu_L1';
    t.P_L2 = Pseu_L2';
    t.P_L5 = Pseu_L5';
%     t.L1 = L1';
%     t.L2 = L2';
%     t.L5 = L5';
    t.C1 = C1';
    t.C2 = C2';
    t.C5 = C5';
    t.P1 = P1';
    t.P2 = P2';
%     t.S1 = S1';
%     t.S2 = S2';
%     t.S5 = S5';
%     t.L1_LLI = L1_LLI';
%     t.L2_LLI = L2_LLI';
%     t.L5_LLI = L5_LLI';
    t.C1_LLI = C1_LLI';
    t.C2_LLI = C2_LLI';
%     t.C5_LLI = C5_LLI';
    t.P1_LLI = P1_LLI';
    t.P2_LLI = P2_LLI';
%     t.S1_LLI = S1_LLI';
%     t.S2_LLI = S2_LLI';
%     t.S5_LLI = S5_LLI';
%     t.L1_str = L1_str';
%     t.L2_str = L2_str';
%     t.L5_str = L5_str';
    t.C1_str = C1_str';
    t.C2_str = C2_str';
%     t.C5_str = C5_str';
    t.P1_str = P1_str';
    t.P2_str = P2_str';
%     t.S1_str = S1_str';
%     t.S2_str = S2_str';
%     t.S5_str = S5_str';

end

function t = read_rinex_obs_legacy( filename)
    t = table();
    fid = fopen(filename);
    if fid == -1
        errordlg(['The file ''' filename ''' does not exist.']);
        return;
    end
    % skip through header
    end_of_header = 0;
    while end_of_header == 0
        current_line = fgetl(fid);
        if contains(current_line,'END OF HEADER')
            end_of_header=1;
        end
    end
    j = 0;
    n = 1;
    while feof(fid) ~= 1
        j = j+1;

        current_line = fgetl(fid);
        
        if (length(current_line) < 33)
            current_line = fgetl(fid);
            current_line = fgetl(fid);
            current_line = fgetl(fid);
            current_line = fgetl(fid);
            current_line = fgetl(fid);
        end
        
        if (length(current_line) < 33)
            F = '?'
        end
        
        Y_t = str2double(current_line(1:3));
        M_t = str2double(current_line(4:6));
        D_t = str2double(current_line(7:9));
        H_t = str2double(current_line(10:12));
        min_t = str2double(current_line(13:15));
        sec_t = str2double(current_line(16:26));
        epoch_flag_t = str2double(current_line(27:29));
        num_sat = str2double(current_line(30:32));
        if(num_sat <=12)
            prns = current_line(33:(32+3*num_sat));
            prn_list = strsplit(prns, 'G');
            prn_array = str2double(prn_list(2:(num_sat+1)));
        else
            prns = current_line(33:(32+3*12));
            prn_list = strsplit(prns, 'G');
            prn_array = str2double(prn_list(2:(12+1)));
            current_line = fgetl(fid);
            prns = current_line(33:(32+3*(num_sat-12)));
            prn_list = strsplit(prns, 'G');
            prn_array = [prn_array str2double(prn_list(2:(num_sat-11)))];  
        end
        
        
        for i = 1:num_sat
            % Broadcast orbit line 1
            current_line = fgetl(fid);
%             if (length(current_line) >= 14) L1(n) = str2double(current_line(1:14)); else L1(n) = NaN; end
%             if (length(current_line) >= 15) L1_LLI(n) = str2double(current_line(15)); else L1_LLI(n) = NaN; end
%             if (length(current_line) >= 16) L1_str(n) = str2double(current_line(16)); else L1_str(n) = NaN; end
%             if (length(current_line) >= 30) L2(n) = str2double(current_line(17:30)); else L2(n) = NaN; end
%             if (length(current_line) >= 31) L2_LLI(n) = str2double(current_line(31)); else L2_LLI(n) = NaN; end
%             if (length(current_line) >= 32) L2_str(n) = str2double(current_line(32)); else L2_str(n) = NaN; end
            if (length(current_line) >= 46) P1(n) = str2double(current_line(33:46)); else P1(n) = NaN; end
            if (length(current_line) >= 47) P1_LLI(n) = str2double(current_line(47)); else P1_LLI(n) = NaN; end
            if (length(current_line) >= 48) P1_str(n) = str2double(current_line(48)); else P1_str(n) = NaN; end
            if (length(current_line) >= 62) P2(n) = str2double(current_line(49:62)); else P2(n) = NaN; end
            if (length(current_line) >= 63) P2_LLI(n) = str2double(current_line(63)); else P2_LLI(n) = NaN; end
            if (length(current_line) >= 64) P2_str(n) = str2double(current_line(64)); else P2_str(n) = NaN; end

            if (length(current_line) >= 78) C1(n) = str2double(current_line(65:78)); else C1(n) = NaN; end
            if (length(current_line) >= 79) C1_LLI(n) = str2double(current_line(79)); else C1_LLI(n) = NaN; end
            if (length(current_line) >= 80) C1_str(n) = str2double(current_line(80)); else C1_str(n) = NaN; end
            
            Y(n) = Y_t;
            M(n) = M_t;
            D(n) = D_t;
            H(n) = H_t;
            min(n) = min_t;
            sec(n) = sec_t;
            epoch_flag(n) = epoch_flag_t;
            prn(n) = prn_array(i);
            
            if (~isnan(C1(n)))
                Pseu_L1(n) = C1(n);
            else
                Pseu_L1(n) = P1(n);
            end
            
            Pseu_L2(n) = P2(n);
            Pseu_L5(n) = NaN;
            
            n = n + 1;
            
        end
    end
    
    t.year = Y';
    t.M = M';
    t.D = D';
    t.H = H';
    t.min = min';
    t.sec = sec';
    [GPSWeek tr] = UT1toGPStime(D,M,Y,H,min,sec);
    t.GPSWeek = GPSWeek;
    t.tr = tr;
    t.ut1time = UT1tosecs(D,M,Y,H,min,sec);
    t.epoch_flag = epoch_flag';
    t.prn = prn';
    t.P_L1 = Pseu_L1';
    t.P_L2 = Pseu_L2';
    t.P_L5 = Pseu_L5';
%     t.L1 = L1';
%     t.L2 = L2';
    t.C1 = C1';
    t.P1 = P1';
    t.P2 = P2';
%     t.L1_LLI = L1_LLI';
%     t.L2_LLI = L2_LLI';
    t.C1_LLI = C1_LLI';
    t.P1_LLI = P1_LLI';
    t.P2_LLI = P2_LLI';
%     t.L1_str = L1_str';
%     t.L2_str = L2_str';
    t.C1_str = C1_str';
    t.P1_str = P1_str';
    t.P2_str = P2_str';

end

function t = read_phone_rinex_obs( filename )
    t = table();
    fid = fopen(filename);
    if fid == -1
        errordlg(['The file ''' filename ''' does not exist.']);
        return;
    end
    % skip through header
    end_of_header = 0;
    while end_of_header == 0
        current_line = fgetl(fid);
        if contains(current_line,'END OF HEADER')
            end_of_header=1;
        end
    end
    j = 0;
    n = 1;
    while feof(fid) ~= 1
        j = j+1;

        current_line = fgetl(fid);
        
%         if (length(current_line) < 35)
%             current_line = fgetl(fid);
%             current_line = fgetl(fid);
%             current_line = fgetl(fid);
%             current_line = fgetl(fid);
%             current_line = fgetl(fid);
%         end
        
        if (length(current_line) < 35)
            F = '?'
        end
        
        Y_t = str2double(current_line(5:6));
        M_t = str2double(current_line(7:9));
        D_t = str2double(current_line(10:12));
        H_t = str2double(current_line(13:15));
        min_t = str2double(current_line(16:18));
        sec_t = str2double(current_line(19:29));
        epoch_flag_t = str2double(current_line(30:32));
        num_sat = str2double(current_line(33:35));  
        
        for i = 1:num_sat
            % Broadcast orbit line
            current_line = fgetl(fid);
            if(current_line(1) == 'G')
                if (length(current_line) >= 3) prn(n) = str2double(current_line(2:3)); else prn(n) = NaN; end
                if (length(current_line) >= 17) C1(n) = str2double(current_line(4:17)); else C1(n) = NaN; end
                if (length(current_line) >= 81) C5(n) = str2double(current_line(68:81)); else C5(n) = NaN; end

                C2(n) = NaN;

                Y(n) = Y_t;
                M(n) = M_t;
                D(n) = D_t;
                H(n) = H_t;
                min(n) = min_t;
                sec(n) = sec_t;
                epoch_flag(n) = epoch_flag_t;

                n = n + 1;
            end
            
        end
    end
    
    t.year = Y';
    t.M = M';
    t.D = D';
    t.H = H';
    t.min = min';
    t.sec = sec';
    [GPSWeek tr] = UT1toGPStime(D,M,Y,H,min,sec);
    t.GPSWeek = GPSWeek;
    t.tr = tr;
    t.ut1time = UT1tosecs(D,M,Y,H,min,sec);
    t.epoch_flag = epoch_flag';
    t.prn = prn';
    t.P_L1 = C1';
    t.P_L2 = C2';
    t.P_L5 = C5';

end
function toc = UT1tosecs(D,M,Y,H,min,sec)

    for i = 1:length(D)
        if M(i) <= 2
            y = Y(i) - 1;
            m = M(i) + 12;
        else
            y = Y(i);
            m = M(i);
        end

        if Y(i) < 1582 || (Y(i) == 1582 && (M(i) < 10 || (M(i)==10 && D(i) <= 4)))
            B = -2 + ((y+4716)/4) - 1179;
        else
            B = y/400 - y/100 + y/4;
        end

        toc(i) = (((365*y - 679004 +floor(B) + floor(30.6001*(m+1)) + D(i))*24 + H(i))*60 + min(i))*60 + sec(i);
    end
    
    toc = toc';
end
function [GPSWeek toc] = UT1toGPStime(D,M,Y,H,min,sec)

    m_to_day = [0,31,31+28,31+28+31,31+28+31+30,31+28+31+30+31,31+28+31+30+31+30,31+28+31+30+31+30+31,31+28+31+30+31+30+31+31,31+28+31+30+31+30+31+31+30,31+28+31+30+31+30+31+31+30+31,31+28+31+30+31+30+31+31+30+31+30];

    y = Y+20;
    d = D - 6;
    m = m_to_day(M);


    B = y/400 - y/100 + y/4;
    
    days = 365*y + floor(B) + m + d;
    GPSWeek = floor(days/7);


    toc = (((365*y + floor(B) + m + d - 7*GPSWeek)*24 + H)*60 + min)*60 + sec;
    
    GPSWeek = GPSWeek';
    toc = toc';
end
