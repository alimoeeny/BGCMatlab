function WindowMotion(src, data, type)
    persistent tslast;
    
    
    if ~isappdata(src, 'MouseState')
        return;
    end
    X = gui.MouseEvent(src, data);

    if strcmp(X.SelectionType,'alt') && X.down
        fprintf('Alt Move (down %d)\n',X.down);
    end
    if X.down == 0
        return;
    end
    ts = now;
    if isempty(tslast) 
        tslast = 0;
    end
    if isfield(X,'count') && X.count > 0
        [swipe, details] = gui.MouseSwipe(X);
        if details.rate < 10
            DATA = GetDataFromFig(src);
            [e, p] = PC.Mouse2Expt(src, data);
            if details.rate < 10 && ishandle(X.h)
                set(X.h,'color','w');
            end
            if (e ~= X.eid || p ~= X.p) && details.meanrate < 0.5
                if DATA.profiling > 1
                    fprintf('New square %.0f,%.0f rate%.2f type %s\n',e,p,details.meanrate,X.SelectionType);
                end
                if ishandle(X.h)
                    PC.DrawBox(e,p, 'celllist',X.h);
                elseif strcmp(X.SelectionType,'open')
                    X.h = PC.DrawBox(e,p, 'celllist', 'color','g');
                else
                    X.h = PC.DrawBox(e,p, 'celllist', 'color','w');
                end
                if ishandle(X.txt)
                    set(X.txt,'string',sprintf('E%dP%d',e,p),'position',[p e-1 0]);
                else
                    X.txt = text(p,e-1,sprintf('E%dP%d',e,p),'color','w');
                end
            end
        end
        t = ts-X.times(1);
        if X.count > 2 && details.rate > 10
            fprintf('swipe %d %.2f %.2f %s rate %.2f\n',X.count,details.dur,details.dx,datestr(X.times(1),'hhmm:ss'),details.dx);
        end
        setappdata(gcbf,'MouseState',X);
    end

tslast = ts;
