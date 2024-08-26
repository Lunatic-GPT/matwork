function anonymize_UID(dname,keep_name)
%anonymize dicom files in dname
dir_str=dir2(dname);


cur_d=cd(dname);
try
    for i=1:length(dir_str)
        
        in=dicominfo(dir_str(i).name);
        
        in=change_fields(in,keep_name);
        x=dicomread(dir_str(i).name);
        dicomwrite(x, dir_str(i).name,in);
    end
catch
    fprintf('error anonymize %s',dir_str(i).name);
end
cd(cur_d);

function in=change_fields(in,keep_name)

fieldName={'StudyDate','SeriesDate','AcquisitionDate','ContentDate','OverlayDate',...
           'CurveDate','AcquisitionDatetime','SeriesTime','AcquisitionTime','ContentTime',...
           'OverlayTime','CurveTime','AccessionNumber','InstitutionName','InstitutionAddress',...
           'ReferringPhysiciansName','ReferringPhysiciansAddress','ReferringPhysiciansTelephoneNumber','ReferringPhysicianIDSequence','InstitutionalDepartmentName',...
           'PhysicianOfRecord','PhysicianOfRecordIDSequence','PerformingPhysiciansName','PerformingPhysicianIDSequence','NameOfPhysicianReadingStudy',...
           'PhysicianReadingStudyIDSequence','IssuerOfPatientID','PatientBirthDate','PatientBirthTime','PatientSex',...
           'OtherPatientIDs','OtherPatientNames','PatientBirthName','PatientAddress','PatientsMothersBirthName',...
           'CountryOfResidence','RegionOfResidence','PatientTelephoneNumbers','StudyID','CurrentPatientLocation',...
           'PatientInstitutionResidence','DateTime','Date','Time','PersonName','PatientAge'...
           'DeviceSerialNumber','PatientWeight','PatientID','PatientSize',...
           'FrameOfReferenceUID','ImplementationClassUID','MediaStorageSOPClassUID',...
           'MediaStorageSOPInstanceUID','SeriesInstanceUID','SOPClassUID','SOPInstanceUID','StudyInstanceUID','TransferSyntaxUID'};
  if ~keep_name
      fieldName{end+1}='PatientName';
  end
       
for  i=1:length(fieldName)
   if isfield(in,fieldName{i})
      fprintf('%s - %s\n',fieldName{i},class(getfield(in,fieldName{i})));
       
      val=getfield(in,fieldName{i});
      if strmatch(class(val),'double')
          in=setfield(in,fieldName{i},99);
      else
          in=setfield(in,fieldName{i},repmat('x',[1,length(val)]));
      end
   end
end

       
       
       