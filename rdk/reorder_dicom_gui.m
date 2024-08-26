function reorder_dicom_gui

for i=1:8
    dcmOrder(sprintf('rdk%d',i),'*.ima',false);
end
