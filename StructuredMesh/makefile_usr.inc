SRC_LOC := ./src
KTH := $(SRC_LOC)/KTH_Toolbox-master
$(OBJDIR)/frame.o             :$(KTH)/driver/frame/frame.f;                         $(F77) -c $(FL2) -I./src $< -o $@
$(OBJDIR)/mntrlog_block.o   :$(KTH)/driver/frame/monitor/mntrlog_block.f MNTRLOGD;               $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrlog.o         :$(KTH)/driver/frame/monitor/mntrlog.f MNTRLOGD;                     $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrtmr_block.o   :$(KTH)/driver/frame/monitor/mntrtmr_block.f MNTRLOGD MNTRTMRD;      $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrtmr.o         :$(KTH)/driver/frame/monitor/mntrtmr.f MNTRLOGD MNTRTMRD FRAMELP;    $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/rprm_block.o      :$(KTH)/driver/frame/runparam/rprm_block.f RPRMD;                     $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/rprm.o            :$(KTH)/driver/frame/runparam/rprm.f RPRMD FRAMELP;                   $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_tools_block.o  :$(KTH)/io/io_tools/io_tools_block.f IOTOOLD;               $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_tools.o        :$(KTH)/io/io_tools/io_tools.f IOTOOLD;                     $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/chkpoint.o        :$(KTH)/io/checkpoint/chkpoint.f CHKPOINTD;                   $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/chkpt_mstp.o      :$(KTH)/io/checkpoint/mstep/chkpt_mstp.f CHKPTMSTPD CHKPOINTD;      $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/conht_tools.o      :$(KTH)/utility/conht/conht_tools.f CONHTD;      $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/IBM_tools.o      :$(SRC_LOC)/IBM_tools.f;      $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/parallfil.o      :$(SRC_LOC)/parallfil.f;      $(F77) -c $(FL2) -I./ $< -o $@

$(OBJDIR)/map2D.o           :$(SRC_LOC)/Adam/map2D.f MAP2D FRAMELP;                  $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/stat.o            :$(SRC_LOC)/Adam/stat.f STATD MAP2D FRAMELP;             $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/stat_IO.o         :$(SRC_LOC)/Adam/stat_IO.f STATD MAP2D FRAMELP;          $(F77) -c $(FL2) -I./ $< -o $@
