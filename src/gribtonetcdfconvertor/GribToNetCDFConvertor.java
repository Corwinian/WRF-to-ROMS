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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;


/**
 *
 * @author corwin
 */
public class GribToNetCDFConvertor 
{
	private static String latName = "lat";
	private static String lonName = "lon";
	
	private static int[] neededValues = {
		1,//Surface pressure
		7,//Terrain height
		13,//Skin potential temperature
		51,//Skin specific humidity
		52,//Skin Relative humidity
		11,//Skin temperature
		33,//10 M u component wind
		34,//10 M v component wind
		204,//Incoming surface shortwave radiation — time-averaged
		205,//Incoming surface longwave radiation - time-averaged
		211,//Outgoing surface shortwave radiation - time-averaged			
		212,//Outgoing surface longwave radiation – time-averaged
		124,//Surface u wind stress
		125,//Surface v wind stress
		122,//Surface sensible heat flux — time-averaged
		155,//Ground heat flux — time-averaged
		121,//Surface latent heat flux — time-averaged
		172,//Surface momentum flux — time-averaged
		91,//Sea ice mask
		92,//Ice thickness
		81,//Land sea mask (land=1, sea=0)
		154,//Accumulated land surface model precipitation
		57,//Accumulated surface evaporation
	};
	
	private static Map<Integer, String> getVariablesByNums(NetcdfFile cdf) throws IOException
	{
		List<Variable> variables = cdf.getVariables();
		
		Map<Integer, String> map = new HashMap<Integer, String>();
		/*я не нашел как сделать через итераторы, поэтому сделал так*/
		for(int i = 0; i < variables.size(); ++i)
		{ 
			Attribute atr = variables.get(i).findAttribute("GRIB_param_number");
			
			if (atr != null)
				map.put((int)atr.getValue(0), variables.get(i).getName());
		}
		
		return map;
	}
	
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
			
			Map<Integer, String> variables = getVariablesByNums(cdf);
			
			double lat[]=loadCoords(cdf, latName);
			double lon[]=loadCoords(cdf, lonName);
			
			createNetCdfFile(outFile, lat, lon, variables);
			
			for(int i =0; i < neededValues.length; ++i)
			{	
				String field = variables.get(neededValues[i]);
				
				if (field == null)
				{
					System.out.println(String.format("Can't find variavle with number: %d", neededValues[i]));
					continue;
				}
				
				System.out.println(field);

				double [][] data = get2DField(cdf, field, lat, lon);

				NetCDFOperator.writeFieldToNetCDF(outFile, field, data);
			}
		}
		finally
		{
			cdf.close();
		}
	}
	
	public static void createNetCdfFile(String dstFile, double latitude[], double longitude[], Map<Integer, String> variables) 
    {
        NetcdfFileWriteable cdf=null;
        int LONNUM = longitude.length;
        int LATNUM = latitude.length;
        
		try
        {
			File outFile = new File(dstFile);
            if (outFile.exists()) 
				outFile.delete();
			
            cdf = NetcdfFileWriteable.createNew(dstFile);
            Dimension lonDim=cdf.addDimension(lonName, LONNUM),
					latDim=cdf.addDimension(latName, LATNUM);
			
            Dimension[] lodim=new Dimension[1];
            lodim[0]=lonDim;
            cdf.addVariable(lonName, DataType.DOUBLE, lodim);
			
            Dimension[] ladim=new Dimension[1];
            ladim[0]=latDim;
            cdf.addVariable(latName, DataType.DOUBLE, ladim);
			
			Attribute att = new Attribute("_FillValue", -9999);
            Attribute att1 = new Attribute("missing_value", -9999);
			
			//TODO: заменить потом на итераторы
			for (int i = 0; i < neededValues.length; ++i)
			{
				Dimension[] todim=new Dimension[2];
				todim[0] = latDim;
				todim[1] = lonDim;
				
				String field = variables.get(neededValues[i]);
				
				if (field == null)
					continue;
				
				Variable var = cdf.addVariable(field, DataType.DOUBLE, todim);
				
				var.addAttribute(att);
				var.addAttribute(att1);
			}
			
			cdf.create();
            cdf.close();
            
			NetCDFOperator.writeTimeToNetCDF(dstFile, lonName, longitude);
            NetCDFOperator.writeTimeToNetCDF(dstFile, latName, latitude);
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
		//String fileIn = "/home/corwin/Dropbox/Учеба/Курсовик/wrfprs.000.grb";
		String fileIn ="/home/corwin/Dropbox/Учеба/Курсовик/wrf/wrfprs_for_roms.003.grb";
		String fileOut = "/home/corwin/Dropbox/Учеба/Курсовик/wrf/wrfprs_for_roms.003.nc";
	   
		File fIn = new File(fileIn);
		//File fOut = new File(fileOut);
		
		
		//GribToNetCDFExtractor.printMetaData(fIn);
		GribToNetCDFConvert(fIn, fileOut);
	}
}
