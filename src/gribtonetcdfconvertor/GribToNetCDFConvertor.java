/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package gribtonetcdfconvertor;

import java.io.*;
import LRFD.Common.*;

import LRFD.DataStructures.*;
import LRFD.Methods.Interpolation;
import LRFD.Db.NetCDFFile.*;

import ucar.nc2.*;
import ucar.ma2.*;
import java.util.List;

/**
 *
 * @author corwin
 */
public class GribToNetCDFConvertor 
{
	private static String latName = "lat";
	private static String lonName = "lon";
	
	private static double[] loadCoords(NetcdfFile cdf, String coordinate) throws IOException
   {
		Variable latitude = cdf.findVariable(coordinate);
		ucar.ma2.Array lats = latitude.read();
		int n = (int)lats.getSize();
		double lat[] = new double[n];
		
		for (int i = 0; i < n; i++)
		{
			lat[i]=lats.getDouble(i);
		}
		return lat;
	}
	

	
	public static double[][] get2DField(NetcdfFile cdf, String variableName, double[] lat, double[] lon) throws IOException
	{
		Variable var=cdf.findVariable(variableName);
		List<Dimension> dims=var.getDimensions();
		int sz=dims.size();
	  
		ucar.ma2.Array varar=var.read();
		Index ind=varar.getIndex();
		if (sz == 2) 
			ind.set(0, 0);
		else if(sz == 3) 
			ind.set(0, 0, 0);
		else if (sz == 4) 
			ind.set(0, 0, 0, 0);
		double res[][]=new double[lat.length][lon.length];
		
		int jf = (sz==2) ? dims.get(0).getLength() : (sz==3 ? dims.get(1).getLength() : dims.get(2).getLength());
		int kf = (sz==2) ? dims.get(1).getLength() : (sz==3 ? dims.get(2).getLength() : dims.get(3).getLength());
		
		for (int j=0; j < jf; j++)
		{
			for (int k=0; k < kf; k++)
			{
				res[j][k]=varar.getDouble(ind);
				ind.incr();//set(i,j,k);
				//System.out.print(String.format("%.2f ", SSTar.getDouble(ind)));
			}
		}
		return res;
	}
	 
	public static void GribToNetCDFConvert(File fileIn, String outFile) throws IOException//, String dstFileName)
	{   
		NetcdfFile cdf = null;
		try
		{
			cdf = NetcdfFile.open(fileIn.getAbsolutePath());
			
			double lat[]=loadCoords(cdf, latName);
			double lon[]=loadCoords(cdf, lonName);
			
			{
				String field = "v_wind";
				System.out.println(field);

				double [][] data = get2DField(cdf, field, lat, lon);

				Write2DFieldToNetcdf(outFile, lat, lon, data, lonName, latName, field);
			}
		}
		finally
		{
			cdf.close();
		}
	}
	
	public static void Write2DFieldToNetcdf(String dstFile, double latitude[], double longitude[], double vals[][],
			String lon_name, String lat_name, String val_name)
	{
		NetcdfFileWriteable cdf=null;
		int LONNUM = longitude.length,
			LATNUM = latitude.length;
		
		try
		{
			cdf = NetcdfFileWriteable.createNew(dstFile);
			Dimension lonDim = cdf.addDimension(lon_name, LONNUM),
					 latDim = cdf.addDimension(lat_name, LATNUM);
			
			Dimension[] lodim = new Dimension[1];
			lodim[0] = lonDim;
			cdf.addVariable(lon_name, DataType.DOUBLE, lodim);
			
			Dimension[] ladim = new Dimension[1];
			ladim[0] = latDim;
			cdf.addVariable(lat_name, DataType.DOUBLE, ladim);
			
			Dimension[] todim = new Dimension[2];
			todim[0] = latDim;
			todim[1] = lonDim;
			Variable var = cdf.addVariable(val_name, DataType.DOUBLE, todim);
            Attribute att = new Attribute("_FillValue", -9999);
            Attribute att1 = new Attribute("missing_value", -9999);
            var.addAttribute(att);
            var.addAttribute(att1);
			
			cdf.create();
			cdf.close();
			
			NetCDFOperator.writeTimeToNetCDF(dstFile, lon_name, longitude);
			NetCDFOperator.writeTimeToNetCDF(dstFile, lat_name, latitude);
			NetCDFOperator.writeFieldToNetCDF(dstFile, val_name, vals);
		}
		catch (Exception er)
		{
			er.printStackTrace();
		}
	}
   
   /**
	* @param args the command line arguments
	*/
	public static void main(String[] args) throws IOException
	{
		String fileIn = "/home/corwin/Dropbox/Учеба/Курсовик/wrf/wrfprs_for_roms.003.grb";
		String fileOut = "/home/corwin/Dropbox/Учеба/Курсовик/wrf/wrfprs_for_roms.003.nc";
	   
		File fIn = new File(fileIn);
		File fOut = new File(fileOut);
		
		if (fOut.exists()) 
			fOut.delete();
		
		//GribToNetCDFExtractor.printMetaData(fIn);
		GribToNetCDFConvert(fIn, fileOut);
	}
}
