close all 
clear all

%% injections
RA_inj = 146.182916666667*pi/180; %  rads  % 9.745527777778  hours
dec_inj = 9.610000000000 *pi/180; % rads
dis_inj = 9400.000000000000 /1000; %MPc
%% read samples 
test_bin_size=0;
plot_key=0;
if (test_bin_size)
data_name=8;
else
    data_name=0:9;
end
for idata=1:length(data_name)    
sample_data=['../../newsamples/posterior_galaxycat_event-' num2str(data_name(idata)) '.dat'];

load_posterior(sample_data);
pos=data;
n_post_sample=length(pos(:,2));
%% read galaxy catalog 
run read_ucng_noplot
ra=ones(length_d,1);
dec=ones(length_d,1);
ra=RAh+RAm./60+RAs./60./60;  %%%% in hour 
north_de=strmatch('+',f_DEd);
south_de=strmatch('-',f_DEd);
dec(north_de)=DEd(north_de)+DEm(north_de)./60+DEm(north_de)./60./60;
dec(south_de)=-DEd(south_de)-DEm(south_de)./60-DEm(south_de)./60./60; %%% in degree
raa=ra.*15.*pi./180; % rads
decc=dec.*pi./180;   % rads
%% 
likelihood_B=abs(M_b); % f function 
%% 3d bins by galaxie locations 
if (test_bin_size==1)
dis_size=0.01:0.01:2;  %%% real size is dis_size*dis
decc_err=0.01:0.01:0.1 ;
raa_err=decc_err ;
else 
  dis_size=0.1;  %%% real size is dis_size*dis
decc_err=0.01;
raa_err=decc_err ;  
end 

for ds=1:length(dis_size)
    for dra=1:length(raa_err)
    
      for gg=1:length_d
    dis_err=dis_size(ds)*dis(gg);
    dis_low=dis(gg)-dis_err;
    dis_high=dis(gg)+dis_err;
    
 raa_low=raa(gg)-raa_err(dra);
 raa_high=raa(gg)+raa_err(dra);
 decc_low=decc(gg)-decc_err(dra);
 decc_high=decc(gg)+decc_err(dra);
ga_pos=find( (dis_low<=pos(:,13)) & (pos(:,13)<=dis_high) & (raa_low<=pos(:,2)) & (pos(:,2)<=raa_high) ...
            & (decc_low<= pos(:,18)) & (pos(:,18)<decc_high)); % wich galaxy is in each cube 
        n_pos_in_bin(gg)=length(ga_pos)/n_post_sample; %% normalized 
         ga_count_g_b(gg)=   n_pos_in_bin(gg)  .*likelihood_B(gg) ./ dis_err;
      % ga_count_g_b(gg)=   n_pos_in_bin(gg)  .*likelihood_B(gg)   
      end    
%post_filename=['post_galaxy_bin.txt'] 
%save(post_filename,  'ga_count_g_b',   '-ascii')
[rank_count_b, id_cont_b]=sort(ga_count_g_b,'descend');
[rank_b, id_b]=sort(likelihood_B,'descend');
n_rank_count_b=length(rank_count_b);
post_total=rank_count_b(1:n_rank_count_b)./(sum(rank_count_b));

length_post_total=length(find(post_total~=0));
post_total_no_zero=post_total(1:length_post_total)
post_no_zero_filename=['post_no_zero_galaxy_bin_.txt'] ;
 ga_count_b_inj_id=find((abs(dis-dis_inj)<0.00001) & (abs(raa-RA_inj)<0.00001) & (abs(decc-dec_inj)<0.0001));

 %%
 if(test_bin_size==0&& plot_key==1)
figure(idata)
plot3(dis_inj,RA_inj,dec_inj,'r*','markersize',6)
hold on
%n_post_use=256; %%% the min number of cout
   plot3(pos(:,13),pos(:,2),pos(:,18),'b.','markersize',4)
  for ff=1:3
 plot3(dis(id_cont_b(ff)),ra(id_cont_b(ff)).*15.*pi./180,dec(id_cont_b(ff)).*pi./180,'rx','markersize',6)
hold on
	 plot3(dis(id_cont_b(ff)),ra(id_cont_b(ff)).*15.*pi./180,dec(id_cont_b(ff)).*pi./180,'ko','markersize',(10-ff)^2)
 end 
plot3(dis_inj,RA_inj,dec_inj,'r*','markersize',6)
grid
xlabel('dis')
ylabel('R.A.')
zlabel('Dec')
   legend('injection','post samples','galaxy in bins')
   post_name=['galaxy_count_b_in_g_bins_event_' num2str(data_name(idata)) '.eps'];
 saveas(gcf,post_name,'epsc2')
 
 end %if test size
 
 post_inj(ds,dra)=ga_count_g_b(ga_count_b_inj_id)./(sum(rank_count_b)) 
 if(length(find(post_total==post_inj(ds,dra)))==1);
 ga_rank(ds,dra)=find(post_total==post_inj(ds,dra));
 else
  ga_rank(ds,dra)=0
 end % end if
 
    end 
end 
end  
%%
if (test_bin_size)
figure(22)
hold on
  plot(dis_size,post_inj(:,1),'r','markersize',6) 
      plot(dis_size,post_inj(:,2),'k','markersize',6) 
       plot(dis_size,post_inj(:,3),'g','markersize',6) 
      plot(dis_size,post_inj(:,5),'b','markersize',6) 
      plot(dis_size,post_inj(:,10),'r:','markersize',6) 
 xlabel('distant size fraction')
ylabel('posterior probability ')
legend('ra err 0.01','ra err 0.02','ra err 0.03','ra err 0.05','ra err 1')
bins_name=['ra_dis_size_pro' num2str(data_name(idata)) '.eps' ]
 saveas(gcf,bins_name,'epsc2')
else
end
