function CheckClusterValues(DATA, C)    if DATA.check.dropi(1) > 0         if C.dropi(3) < DATA.check.dropi(1) && -C.fitdprime(1) > DATA.check.dropi(2)            msg = sprintf('Cluster %d (GM%.2f) May be dropping spikes (dropi %.2f)',C.cluster,C.fitdprime(1),C.dropi(3));            if DATA.interactive > 0                errordlg(msg,'Cluster Warning','modal');            end            PrintMsg(DATA.logname,msg);        end    end        