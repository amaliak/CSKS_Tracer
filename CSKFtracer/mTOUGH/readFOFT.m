function[Observation]=readFOFT()
        obs_count = 0;
        Observation = [];
        for i = 0:9
            filename = ['FOFT_P.0000',num2str(i)];
            if exist(filename,'file')== 2
                [tmp,nFOFT] = FOFT(filename);
                Observation = [Observation;tmp];
                obs_count = obs_count + nFOFT;
            end
        end
        
        for i = 10:99
            filename = ['FOFT_P.000',num2str(i)];
            if exist(filename,'file') == 2
                [tmp,nFOFT] = FOFT(filename);
                Observation = [Observation;tmp(1) tmp(2) tmp(3)];
                obs_count = obs_count + nFOFT;
            end
            if obs_count >= 13
                return;
            end
        end

        for i = 100:999
            filename = ['FOFT_P.00',num2str(i)];
            if exist(filename,'file') == 2
                [tmp,nFOFT] = FOFT(filename);
                Observation = [Observation;tmp(1) tmp(2) tmp(3)];
                obs_count = obs_count + nFOFT;
            end
            if obs_count >= 13
                return;
            end
        end
    
        function [tmp,noFOFT] = FOFT(filename)
            fid = fopen(filename);
            fgetl(fid);
            fgetl(fid);
            tmp = [];
            % read FOFT block for each element
            noFOFT = 0; noGOFT = 0;
            while noFOFT<10 %first key word is FOFT
                s = fgetl(fid); % read line of the time block
                x =sscanf(s,'%*5c %*f %*c %f %*c %f %*c %*f %*c %*f %*c %*f %*c %*f %*c',[1,2]);
                if strcmp(s(1:4),'FOFT')%first key word is FOFT
                    noFOFT = noFOFT + 1;
                    tmp = [tmp;x(1) x(2) 0];
                else % first key word is GOFT
                    noGOFT = noGOFT + 1;
                    tmp(noGOFT,3) = x(2);
                    if noGOFT == noFOFT
                        break;
                    end
                end
            end
        end
    end