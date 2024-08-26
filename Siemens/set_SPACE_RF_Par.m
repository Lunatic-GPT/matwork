function par=set_SPACE_RF_Par(turbo,etl,TE)
par.necho=turbo;
par.esp=etl/turbo;
par.necho_const=round(TE/etl*turbo);
