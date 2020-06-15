% Convergence of IC-GN and analysis of conditioning number of Hessian matrix

ConvItOrNot = USubpb1;
    ConvItOrNot(1:3:end)=ConvItPerEle(:,1);
    ConvItOrNot(2:3:end)=ConvItPerEle(:,1);
    ConvItOrNot(3:3:end)=ConvItPerEle(:,1);
    Plotdisp_show3(ConvItOrNot,coordinatesFEM,elementsFEM)
    figure, imagesc3D(reshape(ConvItPerEle(:,1),MNL));
    
    % %
    CondQSub3by3 = zeros(size(HtempPar,1),1); CondQ = CondQSub3by3;
    HPar = cell(size(HtempPar,2),1); for tempj = 1:size(HtempPar,2), HPar{tempj} = HtempPar(:,tempj); end
    for tempj = 1:size(HtempPar,1)
        waitbar(tempj/size(HtempPar,1))
        HLocal = zeros(12,12);
            HLocal(1:12) = [HPar{1}(tempj); HPar{2}(tempj); HPar{3}(tempj); HPar{4}(tempj); HPar{5}(tempj); HPar{6}(tempj); ...
                HPar{7}(tempj); HPar{8}(tempj); HPar{9}(tempj); HPar{10}(tempj); HPar{11}(tempj); HPar{12}(tempj)];
            HLocal(12+2:12*2) = [HPar{13}(tempj); HPar{14}(tempj); HPar{15}(tempj); HPar{16}(tempj); HPar{17}(tempj); ...
                HPar{18}(tempj); HPar{19}(tempj); HPar{20}(tempj); HPar{21}(tempj); HPar{22}(tempj); HPar{23}(tempj)];
            HLocal(12*2+3:12*3) = [HPar{24}(tempj); HPar{25}(tempj); HPar{26}(tempj); HPar{27}(tempj); HPar{28}(tempj); ...
                HPar{29}(tempj); HPar{30}(tempj); HPar{31}(tempj); HPar{32}(tempj); HPar{33}(tempj)];
            HLocal(12*3+4:12*4) = [HPar{34}(tempj); HPar{35}(tempj); HPar{36}(tempj); HPar{37}(tempj); HPar{38}(tempj); ...
                HPar{39}(tempj); HPar{40}(tempj); HPar{41}(tempj); HPar{42}(tempj)];
            HLocal(12*4+5:12*5) = [HPar{43}(tempj); HPar{44}(tempj); HPar{45}(tempj); HPar{46}(tempj); HPar{47}(tempj); ...
                HPar{48}(tempj); HPar{49}(tempj); HPar{50}(tempj)];
            HLocal(12*5+6:12*6) = [HPar{51}(tempj); HPar{52}(tempj); HPar{53}(tempj); HPar{54}(tempj); HPar{55}(tempj); ...
                HPar{56}(tempj); HPar{57}(tempj)];
            HLocal(12*6+7:12*7) = [HPar{58}(tempj); HPar{59}(tempj); HPar{60}(tempj); HPar{61}(tempj); HPar{62}(tempj); ...
                HPar{63}(tempj)];
            HLocal(12*7+8:12*8) = [HPar{64}(tempj); HPar{65}(tempj); HPar{66}(tempj); HPar{67}(tempj); HPar{68}(tempj)];
            HLocal(12*8+9:12*9) = [HPar{69}(tempj); HPar{70}(tempj); HPar{71}(tempj); HPar{72}(tempj)];
            HLocal(12*9+10:12*10) = [HPar{73}(tempj); HPar{74}(tempj); HPar{75}(tempj)];
            HLocal(12*10+11:12*11) = [HPar{76}(tempj); HPar{77}(tempj)];
            HLocal(12*12) = [HPar{78}(tempj)];
        
        HLocal = HLocal'+HLocal-diag(diag(HLocal));
        
        HLocalSub = zeros(3,3);
        HLocalSub(1:3) = [HPar{1}(tempj); HPar{2}(tempj); HPar{3}(tempj)];
        HLocalSub(5:6) = [HPar{5}(tempj); HPar{6}(tempj)]; 
        HLocalSub(9)   = [HPar{7}(tempj)];
        HLocalSub = HLocalSub'+HLocalSub-diag(diag(HLocalSub));
        
        try
            CondQ(tempj) = cond(inv(HLocal));     
            CondQSub3by3(tempj) = cond(inv(HLocalSub));    
        catch
           CondQ(tempj) = NaN; 
        end
    end
     figure, imagesc3D(reshape(log(CondQ),MNL)); caxis([5,8]);
     figure, imagesc3D(reshape(log(CondQSub3by3),MNL)); caxis([5,8]);