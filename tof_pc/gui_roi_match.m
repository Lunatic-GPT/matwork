function p= gui_roi_match()

p = parameter('gui to match clusters');

p=add(p,'filename','mask 1','*MID??_fl*\interp*\*vmask.mat');
p=add(p,'filename','mask 2','*MID??_fl*\interp*\*vmask.mat');

p=add(p,'filename','underlay 1','');
p=add(p,'filename','underlay 2','');


p=add(p,'float','max dist (mm)',0.3);
p=add(p,'float','vox size (mm)',0.1563);

p=add(p,'string','suffix','_match');
p=add(p,'button','match','gui_roi_match_callback(params)');
p = add(p,'button','Close','close');

if nargout==0

parametergui(p);
end

