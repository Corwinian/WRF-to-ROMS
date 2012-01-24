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
	private static String timeName = "time";
	private static String time1Name = "time1";
	private static String levelName = "height_above_ground";
	private static String level1Name = "height_above_ground1";
	
	
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
	
	private static double[] loadCoords(NetcdfFile cdf, String coordinate, double min, double max) throws IOException, Exception
	{
		return loadCoords(cdf, coordinate, min, max, false);
	}
	
	private static double[] loadCoords(NetcdfFile cdf, String coordinate, double min, double max, boolean inv) throws IOException, Exception
	{
		Variable var = cdf.findVariable(coordinate);
		
		Range range=getRange(var, min, max, false);
		return selectRange(var, range, inv);
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
	
	public static Range getRange(Variable Var, double minValue, double maxValue) throws Exception
	{
		return getRange(Var, minValue, maxValue, false);
	}
	
	public static Range getRange(Variable Var, double minValue, double maxValue, boolean wide) throws Exception
	{
		ucar.ma2.Array tvar=Var.read();
		int iMin=0, iMax=0;
		if (Var.getRank()!=1) throw new Exception("Error in dimensions. Rank != 1");

		int n=Var.getDimension(0).getLength();
		double minVal=Double.MAX_VALUE;
		double maxVal=-Double.MAX_VALUE;

		for (int i=0; i < n; i++)
		{
			double timevar = tvar.getDouble(i);
			
			if (timevar >= minValue && timevar <= maxValue)
			{
				if (timevar < minVal)
				{
					iMin=i;
					minVal=timevar;
				}
				if (timevar > maxVal)
				{
					iMax=i;
					maxVal=timevar;
				}
			}
		}
		if (iMax<iMin)
		{
			int temp=iMax;
			iMax=iMin;
			iMin=temp;
		}
		if (minVal!=minValue || maxVal!=maxValue)
		{
			iMin--;
			iMax++;
			if (iMax>=n) 
				iMax=n;
			if (iMin<0) 
				iMin=0;
		}
		if (wide)
		{
			iMin-=3;
			iMax+=3;
			if (iMax>=n) 
				iMax=n;
			if (iMin<0) 
				iMin=0;
		}
		return new Range(iMin, iMax);
	}
	
	public static double[] selectRange(Variable Var, Range range) throws IOException
	{ return selectRange(Var, range, false);}
	
	public static double[] selectRange(Variable Var, Range range, boolean inv) throws IOException
	{
		double res[]=new double[range.last()-range.first()+1];
		ucar.ma2.Array tvar=Var.read();
		int first=range.first();
		
		for (int i=0; i<res.length; i++)
		{
			res[i] = inv ? tvar.getDouble(res.length - 1 - i + first) : tvar.getDouble(i+first);
		}
		
		return res;
	}
	
	public static Data3DField getFieldFromSRCFile(NetcdfFile src, String FieldName, GeoRectangle boundaries, double timeMin, double timeMax) throws Exception
	{
		return getFieldFromSRCFile(src, FieldName, boundaries, timeMin, timeMax, false);
	}
	public static Data3DField getFieldFromSRCFile(NetcdfFile src, String FieldName, GeoRectangle boundaries, double timeMin, double timeMax, boolean wide) throws Exception
	{
		Variable var=src.findVariable(FieldName);
		if (!(var.getRank()==3 || (var.getRank()==4 && var.getDimension(1).getLength()==1))) 
			throw new Exception("Error in Dimensions, Rank of dimensions for variable "+var.getName()+" must be 3.");
		
		String lvName = (FieldName.equals("Specific_humidity") || FieldName.equals("Relative_humidity")) ?
				level1Name : levelName ;
		String tmName = (FieldName.equals("Evaporation") || FieldName.equals("Land_Surface_Precipitation_Accumulation_LSPA")) ?
				time1Name : timeName;
		
		int nv=var.getRank();
		Variable timeVar=src.findVariable(tmName);
		Variable lonVar=src.findVariable(lonName);
		Variable latVar=src.findVariable(latName);

		Range timeRange=getRange(timeVar, timeMin, timeMax);
		double time[]=selectRange(timeVar, timeRange);

		Range lonRange=getRange(lonVar, boundaries.north_west.lon, boundaries.south_east.lon, wide);
		double lon[]=selectRange(lonVar, lonRange);

		Range latRange=getRange(latVar, boundaries.south_east.lat, boundaries.north_west.lat, wide);
		double lat[]=selectRange(latVar, latRange);

		ArrayList<Range> lr=new ArrayList<Range>();
		lr.add(var.findDimensionIndex(tmName), timeRange);
		if (nv==4) 
		{
			Range levelRange=new Range(0, 0);
		
			lr.add(var.findDimensionIndex(lvName), levelRange);
		}
		
		lr.add(var.findDimensionIndex(latName), latRange);
		lr.add(var.findDimensionIndex(lonName), lonRange);

		ucar.ma2.Array tvar=var.read(lr);
		long nt=tvar.getSize();
		if (nt!=lon.length*lat.length*time.length)
		{
			System.out.println("Error in dimensions");
		}

		int[] dims=new int[nv];
		dims[var.findDimensionIndex(tmName)]=time.length;
		if (nv==4) 
			dims[var.findDimensionIndex(lvName)]=1;
		dims[var.findDimensionIndex(latName)]=lat.length;
		dims[var.findDimensionIndex(lonName)]=lon.length;
		
		double res[][][] = (nv != 4) ? new double[dims[0]][dims[1]][dims[2]]:new double[dims[0]][dims[2]][dims[3]];
		double scaleFactor=1;
		double addOffset=0;
		boolean kridControl=false;
		double minV=0;
		double maxV=0;
		double mis=Double.MAX_VALUE;
		double fill=Double.MAX_VALUE;
		try
		{
			Attribute SF=var.findAttribute("scale_factor");
			Attribute AO=var.findAttribute("add_offset");

			scaleFactor=SF.getNumericValue().doubleValue();
			addOffset=AO.getNumericValue().doubleValue();
			Attribute VR=var.findAttribute("valid_range");
			minV=VR.getNumericValue(0).doubleValue();
			maxV=VR.getNumericValue(1).doubleValue();
			minV=minV*scaleFactor+addOffset;
			maxV=maxV*scaleFactor+addOffset;
			kridControl=true;
			Attribute MIS=var.findAttribute("missing_value");
			Attribute FILL=var.findAttribute("_FillValue");
			mis=MIS.getNumericValue().doubleValue();
			fill=FILL.getNumericValue().doubleValue();
			//mis=mis*scaleFactor+addOffset;
			//fill=fill*scaleFactor+addOffset;
		}
		catch (Exception er){}
		
		Index ima=tvar.getIndex();
		int first=dims[0];
		int second = (nv==4) ? dims[2]:dims[1];
		int third = (nv==4) ? dims[3]:dims[2];
		
		for (int i=0; i<first; i++)
		{
			for (int j=0; j<second; j++)
			{
				for (int k=0; k<third; k++)
				{
					if (nv==3) 
						ima.set(i, j, k);
					else if (nv==4) 
						ima.set(i, 0, j, k);
					res[i][j][k]=tvar.getDouble(ima);
					
					if (res[i][j][k]==mis)
					{
						System.out.println("Missing value = "+res[i][j][k]+" in position "+i+"; "+j+"; "+k+"; replaced by "+fill);
						res[i][j][k]=fill;
					}
					
					res[i][j][k]=res[i][j][k]*scaleFactor+addOffset;
					if (kridControl)
					{					  
						if (res[i][j][k]<minV || res[i][j][k]>maxV)
						{
							System.out.println("Illegar value = "+res[i][j][k]+" in position "+i+"; "+j+"; "+k+";");
							//throw new Exception("Illegar value in position "+i+"; "+j+"; "+k+";");
						}
					}
				}
			}
		}
		return new Data3DField(res, lat, lon, time);
	}
	
	public static Data3DField InterpolateField(Data3DField field, double[] res_lat, double[] res_lon) throws Exception
	{
		Data3DField resflx=new Data3DField();
		resflx.lat=res_lat;
		resflx.lon=res_lon;
		resflx.time=field.time;
		resflx=Interpolation.BilinearInterpolation(field, resflx);
		return resflx;
	}
	
	public static void GribToNetCDFConvert(File fileIn, String outFile) throws IOException, Exception //, String dstFileName)
	{   
		NetcdfFile cdf = null;
		try
		{ 
			GeoRectangle gr = new GeoRectangle((float)(43-4), (float)(63+4), (float)(129-4), (float)(160+4));
			cdf = NetcdfFile.open(fileIn.getAbsolutePath());

			Map<Integer, String> variables = getVariablesByNums(cdf);
			
			double lat[]=loadCoords(cdf, latName, gr.south_east.lat+4, gr.north_west.lat-4, true);
			double lon[]=loadCoords(cdf, lonName, gr.north_west.lon+4, gr.south_east.lon-4);
			double time[]=loadCoords(cdf, timeName);
			double time1[]=loadCoords(cdf, time1Name);
			
			createNetCdfFile(outFile, lat, lon, time, variables);

			for(int i =0; i < neededValues.length; ++i)
			{	
				String field = variables.get(neededValues[i]);
				
				if (field == null)
				{
					System.out.println(String.format("Can't find variavle with number: %d", neededValues[i]));
					continue;
				}
				
				System.out.println(field);
				Data3DField SST= getFieldFromSRCFile(cdf, field,gr, time[0], time[time.length -1]);
				SST.InverseLatIfNecessary();
				SST = InterpolateField(SST, lat, lon) ;
				SST.InverseLatIfNecessary();

				NetCDFOperator.writeFieldToNetCDF(outFile, field,  SST.data);
			}
		}
		finally
		{
			cdf.close();
		}
	}
	
	public static void createNetCdfFile(String dstFile, double latitude[], double longitude[],
			double time[], Map<Integer, String> variables) 
	{
		NetcdfFileWriteable cdf=null;
		int LONNUM = longitude.length;
		int LATNUM = latitude.length;
		int TIMENUM = time.length;
		
		try
		{
			File outFile = new File(dstFile);
			if (outFile.exists()) 
				outFile.delete();
			
			cdf = NetcdfFileWriteable.createNew(dstFile);
			Dimension lonDim=cdf.addDimension(lonName, LONNUM),
					latDim=cdf.addDimension(latName, LATNUM),
					timeDim=cdf.addDimension(timeName, TIMENUM);
			
			Dimension[] lodim=new Dimension[1];
			lodim[0]=lonDim;
			cdf.addVariable(lonName, DataType.DOUBLE, lodim);
			
			Dimension[] ladim=new Dimension[1];
			ladim[0]=latDim;
			cdf.addVariable(latName, DataType.DOUBLE, ladim);
			
			Dimension[] tidim=new Dimension[1];
			tidim[0]=timeDim;
			cdf.addVariable(timeName, DataType.DOUBLE, tidim);
			
			Attribute att = new Attribute("_FillValue", -9999);
			Attribute att1 = new Attribute("missing_value", -9999);
			
			//TODO: заменить потом на итераторы
			for (int i = 0; i < neededValues.length; ++i)
			{
				Dimension[] todim=new Dimension[3];
				todim[0] = timeDim;
				todim[1] = latDim;
				todim[2] = lonDim;
				
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
			NetCDFOperator.writeTimeToNetCDF(dstFile, timeName, time);
		}
		catch (Exception er)
		{
			er.printStackTrace();
		}
	}
	
   /**
	* @param args the command line arguments
	*/
	public static void main(String[] args) throws IOException, Exception
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
