function [barcode_syn,barcode_rec,a,m0,community_name]=generate_mix_invasion(k,num_strain,community,invader)

barcode_syn = [community.barcode_syn_record(:,1:num_strain-1,k) invader.barcode_syn_record(:,num_strain,k)];
barcode_rec = [community.barcode_rec_record(:,1:num_strain-1,k) invader.barcode_rec_record(:,num_strain,k)];
a = [community.a_record(k,1:num_strain-1) invader.a_record(k,num_strain)];
m0 = [community.m0_record(k,1:num_strain-1) invader.m0_record(k,num_strain)];
a = [a;1-a];
community_name=[inputname(4) '_invade_' inputname(3)];
end