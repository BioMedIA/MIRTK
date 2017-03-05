/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "mirtk/DataOp.h"

#include "mirtk/Path.h"
#include "mirtk/Memory.h"
#include "mirtk/BaseImage.h"
#include "mirtk/ImageAttributes.h"

#include "mirtk/ImageConfig.h"
#if MIRTK_Image_WITH_VTK
  #include "vtkSmartPointer.h"
  #include "vtkDataSetReader.h"
  #include "vtkDataSetWriter.h"
  #include "vtkXMLGenericDataObjectReader.h"
  #include "vtkXMLDataSetWriter.h"
  #include "vtkImageData.h"
  #include "vtkPointSet.h"
  #include "vtkPointData.h"
  #include "vtkCellData.h"
  #include "vtkFloatArray.h"
  #include "mirtk/Vtk.h"
#endif


namespace mirtk { namespace data {


// -----------------------------------------------------------------------------
DataFileType FileType(const char *name)
{
  string ext = Extension(name);
  if (ext == ".vtk") return LEGACY_VTK;
  if (ext.length() == 4 && ext[0] == '.' && ext[1] == 'v') return XML_VTK;
  return IMAGE;
}

// -----------------------------------------------------------------------------
int Read(const char *name, UniquePtr<double[]> &data, int *dtype, ImageAttributes *attr,
         #if MIRTK_Image_WITH_VTK
           vtkSmartPointer<vtkDataSet> *dataset,
         #else
           void *,
         #endif
         const char *scalars_name, bool cell_data)
{
  int n = 0;
  data.reset();
  if (attr) *attr = ImageAttributes();
#if MIRTK_Image_WITH_VTK
  if (dataset) *dataset = nullptr;
#endif // MIRTK_Image_WITH_VTK
  int ftype = FileType(name);
  switch (ftype) {
#if MIRTK_Image_WITH_VTK
    case LEGACY_VTK:
    case XML_VTK: {
      vtkSmartPointer<vtkDataArray> scalars;
      if (ftype == LEGACY_VTK) {
        vtkSmartPointer<vtkDataSetReader> reader;
        reader = vtkSmartPointer<vtkDataSetReader>::New();
        reader->SetFileName(name);
        reader->Update();
        vtkDataSet * const output = reader->GetOutput();
        if (output) {
          vtkDataSetAttributes *arrays = output->GetPointData();
          if (cell_data) arrays = output->GetCellData();
          if (scalars_name && scalars_name[0] != '\0') {
            scalars = arrays->GetArray(scalars_name);
          } else {
            scalars = arrays->GetScalars();
          }
        }
        if (dataset) *dataset = output;
      } else {
        vtkSmartPointer<vtkXMLGenericDataObjectReader> reader;
        reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
        reader->SetFileName(name);
        reader->Update();
        vtkDataSet *output = vtkDataSet::SafeDownCast(reader->GetOutput());
        if (output) {
          vtkDataSetAttributes *arrays = output->GetPointData();
          if (cell_data) arrays = output->GetCellData();
          if (scalars_name && scalars_name[0] != '\0') {
            scalars = arrays->GetArray(scalars_name);
          } else {
            scalars = arrays->GetScalars();
          }
        }
        if (dataset) *dataset = output;
      }
      if (!scalars) {
        cerr << "Failed to read VTK dataset! Type is either not supported or dataset has no scalar " << (cell_data ? "cell" : "point") << " data." << endl;
        cerr << "Use -point/-cell-data option to specify the name of a point/cell data array to use instead." << endl;
        exit(1);
      }
      if (dtype) *dtype = FromVTKDataType(scalars->GetDataType());
      n = static_cast<int>(scalars->GetNumberOfTuples()) * scalars->GetNumberOfComponents();
      if (n == 0) {
        cerr << "VTK dataset has empty scalar " << (cell_data ? "cell" : "point") << " data!" << endl;
        exit(1);
      }
      data.reset(Allocate<double>(n));
      double *p = data.get();
      for (vtkIdType i = 0; i < scalars->GetNumberOfTuples(); ++i) {
        for (int j = 0; j < scalars->GetNumberOfComponents(); ++j, ++p) {
          *p = scalars->GetComponent(i, j);
        }
      }
    } break;
#else // MIRTK_Image_WITH_VTK
    case LEGACY_VTK:
    case XML_VTK:
      cerr << "Cannot read VTK files when Image module not built WITH_VTK!" << endl;
      exit(1);
#endif // MIRTK_Image_WITH_VTK
    case IMAGE: {
      UniquePtr<BaseImage> image(BaseImage::New(name));
      if (attr) *attr = image->Attributes();
      if (dtype) *dtype = image->GetDataType();
      n = image->NumberOfVoxels();
      data.reset(Allocate<double>(n));
      for (int i = 0; i < n; ++i) {
        data[i] = image->GetAsDouble(i);
      }
    } break;
    default:
      cerr << "Unsupported input data file: " << name << endl;
      exit(1);
  }
  return n;
}

// -----------------------------------------------------------------------------
void Write::Process(int n, double *data, bool *)
{
  if (_FileName.empty()) {
    cerr << "Output file name not set!" << endl;
    exit(1);
  }
  int type = FileType(_FileName.c_str());
  switch (type) {
#if MIRTK_Image_WITH_VTK
    case LEGACY_VTK:
    case XML_VTK: {
      if (!_DataSet) {
        cerr << "Cannot write data sequence to VTK file when input was not a VTK file itself!" << endl;
        exit(1);
      }
      vtkDataArray *input_scalars;
      vtkDataSetAttributes *arrays = _DataSet->GetPointData();
      if (_AsCellData)      arrays = _DataSet->GetCellData();
      if (_ArrayName.empty()) input_scalars = arrays->GetScalars();
      else                    input_scalars = arrays->GetArray(_ArrayName.c_str());
      if (input_scalars == nullptr) {
        cerr << "Invalid output array name " << _ArrayName << endl;
        exit(1);
      }
      vtkIdType ntuples = input_scalars->GetNumberOfTuples();
      int       m       = input_scalars->GetNumberOfComponents();
      if (n != static_cast<int>(ntuples * m)) {
        cerr << "Cannot write data sequence to file! Length of data sequence changed." << endl;
        exit(1);
      }
      vtkSmartPointer<vtkDataArray> output_scalars = NewVtkDataArray(ToVTKDataType(_DataType));
      if (!_OutputName.empty()) {
        output_scalars->SetName(_OutputName.c_str());
      } else if (!_ArrayName.empty()) {
        output_scalars->SetName(_ArrayName.c_str());
      }
      output_scalars->SetNumberOfComponents(m);
      output_scalars->SetNumberOfTuples(ntuples);
      for (int j = 0; j < m; ++j) {
        output_scalars->SetComponentName(j, input_scalars->GetComponentName(j));
      }
      for (vtkIdType i = 0; i < ntuples; ++i) {
        for (int j = 0; j < m; ++j, ++data) {
          output_scalars->SetComponent(i, j, *data);
        }
      }
      if (input_scalars == arrays->GetScalars() && _OutputName.empty()) {
        arrays->SetScalars(output_scalars);
      } else {
        if (!_ArrayName.empty() && (_OutputName.empty() || _OutputName == _ArrayName)) {
          arrays->RemoveArray(_ArrayName.c_str());
        }
        arrays->AddArray(output_scalars);
      }
      if (type == LEGACY_VTK) {
        vtkSmartPointer<vtkDataSetWriter> writer;
        writer = vtkSmartPointer<vtkDataSetWriter>::New();
        writer->SetFileName(_FileName.c_str());
        SetVTKInput(writer, _DataSet);
        writer->Write();
      } else {
        vtkSmartPointer<vtkXMLDataSetWriter> writer;
        writer = vtkSmartPointer<vtkXMLDataSetWriter>::New();
        writer->SetFileName(_FileName.c_str());
        SetVTKInput(writer, _DataSet);
        writer->Write();
      }
    } break;
#else
    case LEGACY_VTK:
    case XML_VTK:
      cerr << "Cannot write data series to VTK file when Image module not built WITH_VTK" << endl;
      exit(1);
#endif // MIRTK_Image_WITH_VTK
    case IMAGE: {
      if (_Attributes.NumberOfLatticePoints() == 0) {
        cerr << "Cannot write data series to image file when input was not an image file itself!" << endl;
        exit(1);
      }
      if (_Attributes.NumberOfLatticePoints() != n) {
        cerr << "Cannot write data series to file! Length of data series changed." << endl;
        exit(1);
      }
      UniquePtr<BaseImage> image(BaseImage::New(_DataType));
      image->Initialize(_Attributes);
      for (int i = 0; i < n; ++i) image->PutAsDouble(i, data[i]);
      image->Write(_FileName.c_str());
    } break;
    default:
      cerr << "Unknown output file type: " << _FileName << endl;
      exit(1);
  }
}


} } // namespace mirtk::data
