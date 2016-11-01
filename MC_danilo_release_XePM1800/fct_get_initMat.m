function [texit v0mat r0mat] = fct_get_initMat(texit,v0x_sample,v0y_sample,v0z_sample,r0x_sample,r0y_sample,r0z_sample)

v0mat = [v0x_sample' v0y_sample' v0z_sample'];
r0mat = [r0x_sample' r0y_sample' r0z_sample'];

end