clear
close all

t = readtable("IF_output.csv");

[n_f,IF] = downsample(1,t,1000);

writetable(IF,"IF_output_ds.csv")



function [n_f,IF] = downsample(n_i,IF,ds)
    
    IFnew = table();
    leng = length(IF.X);
    n = n_i;
    s = 0;
    m = 1;
    
    while (n<leng)
        
        time = IF.time(n);
        t_check = time;
        
        while (t_check == time)&&(n<leng)

            if(mod(s,ds) == 0)
                Svid(m) = IF.Svid(n);
                timeIF(m) = IF.time(n);
                A(m) = IF.A(n);
                P(m) = IF.P(n);
                X(m) = IF.X(n);
                Y(m) = IF.Y(n);
                Z(m) = IF.Z(n);
                x(m) = IF.x(n);
                y(m) = IF.y(n);
                z(m) = IF.z(n);
                lat_bot(m) = IF.lat_bot(n);
                long_bot(m) = IF.long_bot(n);
                lat_top(m) = IF.lat_top(n);
                long_top(m) = IF.long_top(n);
                x_bot_ECEF(m) = IF.x_bot_ECEF(n);
                y_bot_ECEF(m) = IF.y_bot_ECEF(n);
                z_bot_ECEF(m) = IF.z_bot_ECEF(n);
                x_top_ECEF(m) = IF.x_top_ECEF(n);
                y_top_ECEF(m) = IF.y_top_ECEF(n);
                z_top_ECEF(m) = IF.z_top_ECEF(n);

                m = m + 1;
            end

            n = n + 1;
            t_check = IF.time(n);
            
        end
        
        s = s+1;
        
    end
    
    n_f = n;

    if(m>1)
        IFnew.Svid = Svid';
        IFnew.time = timeIF';
        IFnew.A = A';
        IFnew.P = P';
        IFnew.X = X';
        IFnew.Y = Y';
        IFnew.Z = Z';
        IFnew.x = x';
        IFnew.y = y';
        IFnew.z = z';
        IFnew.lat_bot = lat_bot';
        IFnew.long_bot = long_bot';
        IFnew.lat_top = lat_top';
        IFnew.long_top = long_top';
        IFnew.x_bot_ECEF = x_bot_ECEF';
        IFnew.y_bot_ECEF = y_bot_ECEF';
        IFnew.z_bot_ECEF = z_bot_ECEF';
        IFnew.x_top_ECEF = x_top_ECEF';
        IFnew.y_top_ECEF = y_top_ECEF';
        IFnew.z_top_ECEF = z_top_ECEF';

    else
        IFnew.Svid = 0;
        IFnew.time = 0;
        IFnew.A = 0;
        IFnew.P = 0;
        IFnew.X = 0;
        IFnew.Y = 0;
        IFnew.Z = 0;
        IFnew.x = 0;
        IFnew.y = 0;
        IFnew.z = 0;
        IFnew.lat_bot = 0;
        IFnew.long_bot = 0;
        IFnew.lat_top = 0;
        IFnew.long_top = 0;
        IFnew.x_bot_ECEF = 0;
        IFnew.y_bot_ECEF = 0;
        IFnew.z_bot_ECEF = 0;
        IFnew.x_top_ECEF = 0;
        IFnew.y_top_ECEF = 0;
        IFnew.z_top_ECEF = 0;
    end
    
    
%     while (t_check == time)&&(n<leng)
%         X(s) = t.X(n);
%         Y(s) = t.Y(n);
%         Z(s) = t.Z(n);
%         B(s) = t.B(n);
%         
%         f1_2 = (1575.42*(10^6))^2;
%         f2_2 = (1227.60*(10^6))^2;
% 
%         
%         a = f1_2/(f1_2 - f2_2);
%         b = f2_2/(f1_2 - f2_2);
%         c = f1_2*f2_2/(f2_2 - f1_2);
%         AIF(s) = c*(t.P_L1(n) - t.P_L2(n));
%         P(s) = a*t.P_L1(n) - b*t.P_L2(n);
%         timeIF(s) = t.tr(n);
%         s = s+1;
%         svid(s) = t.prn(n);
%         
%         n = n + 1;
%         t_check = t.tr(n);
%     end
%     
%     ID = svid(2:end);
%     n_f = n;
%     
%     IF.Svid = ID';
%     IF.time = timeIF';
%     IF.A = AIF';
%     IF.P = P';
%     IF.X = X';
%     IF.Y = Y';
%     IF.Z = Z';
      
end
