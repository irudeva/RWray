% region = {'region   '}   
for ir = 1:length(lat1(1,:))
        region(ir,:) = {sprintf('%d-%d',lat1(ir),lat2(ir))};
        
        for is=5:16
        
        for ic = 1:length(reg(:,1))
            corr = corrcoef(100*squeeze(Kwhite(ir,is,:)/ Ktotal(ir,is,fyr1-syr+1)),squeeze(treg(ic,is,:)));
%           tbl(ir+1,ic+1) = {corr(1,2)};
            corrtbl(is,ir,ic) = corr(1,2);
        end
        end
        
end

varNames = {'region   '};
    for ic = 1:length(reg(:,1))
        varNames(1,ic+1) = {reg(ic,:)};
    end

 
 for is=5:16
     varNames(1,1) = {ssn(is,:)};
     tmp =  squeeze(corrtbl(is,:,:));
     T = table(region, tmp(:,1));
     for ir =2:length(reg(:,1))
         T(:,ir+1)= table(tmp(:,ir));
     end
     T.Properties.VariableNames = varNames;
     
     writetable(T,'../output/Ks/Ks_Treg_corr.xls','Sheet',1,'Range',sprintf('A%d:I%d',(is-5)*5+1,(is-5)*5+5))
 end

        
