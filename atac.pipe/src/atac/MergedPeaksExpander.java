package atac;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;
import java.util.stream.*;

public class MergedPeaksExpander {
	
	String 	MergedIntervalFilePath;
	String 	SampleNameListString;
	String 	OutFilePath;
	boolean isNarrow;
	int 	minRep;
	
	public static void main(String[] args) { 
			boolean DEBUG = true;
			
			args = new String[10];
			
			if(DEBUG) {
				args[0] = "-m";
				args[1] = "/home/aabdalla/Desktop/NextGen2/scratch/data/projects/kramann/jack.3/4_peaks/withDup/FS3.merged.txt";
				args[2] = "-s";
				args[3] = "CD24-1.3,CD24-2.3,CD24-3.3,CD13-1.3,CD13-2.3,CD13-3.3";
				args[4] = "-o";
				args[5] = "/home/aabdalla/Desktop/NextGen2/scratch/data/projects/kramann/jack.3/4_peaks/withDup/FS3.merged.expanded.txt";
				args[6] = "-n";
				args[7] = "true";
				args[8] = "-r";
				args[9] = "1";
			}
		
			// Mandatory arguments
			String MergedIntervalFilePath = "";
			String SampleNameListString = "";
			String OutFilePath = "";
			// Optional arguments (+ defaults)
			boolean isNarrow = true;
			int minRep = 1;
			
			// Scan arguments.
			for(int i = 0; i < args.length; i+=2) {
				switch(args[i]) {
					case "-m"				: MergedIntervalFilePath=args[i+1]			; break;
					case "--merged-peaks"	: MergedIntervalFilePath=args[i+1]			; break;
					case "-s"				: SampleNameListString=args[i+1]			; break;
					case "--samples"		: SampleNameListString=args[i+1]			; break;
					case "-o"				: OutFilePath=args[i+1]						; break;
					case "--output-file"	: OutFilePath=args[i+1]						; break;
					case "-n"				: isNarrow=Boolean.parseBoolean(args[i+1])	; break;
					case "--isNarrow"		: isNarrow=Boolean.parseBoolean(args[i+1])	; break;
					case "-r"				: minRep=Integer.parseInt(args[i+1])		; break;
					case "--min-rep"		: minRep=Integer.parseInt(args[i+1])		; break;
				}
			}
			
			// Run expansion procedure.
			MergedPeaksExpander mpe = new MergedPeaksExpander(MergedIntervalFilePath, SampleNameListString, OutFilePath, isNarrow, minRep);
			mpe.macs2_merged_expand();
	}
		
	public MergedPeaksExpander(String MergedIntervalFilePath, String SampleNameListString, String OutFilePath, boolean isNarrow, int minRep) {
		this.MergedIntervalFilePath = MergedIntervalFilePath;
		this.SampleNameListString = SampleNameListString;
		this.OutFilePath = OutFilePath;
		this.isNarrow = isNarrow;
		this.minRep = minRep;
	}

		/**
		 *  MergedIntervalFile is file created using commands below:
		 *  sort -k1,1 -k2,2n <MACS_NARROWPEAK_FILE_LIST> | 
		 *  mergeBed -c 2,3,4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse 
		 *  > merged_peaks.txt
		 */
			
	public void macs2_merged_expand() {
			
		File MergedIntervalFile = new File(MergedIntervalFilePath);
		File OutFile = new File(OutFilePath);
			 OutFile.getParentFile().mkdirs();
		int  totalOutIntervals = 0;
		List<String> SampleNameList = Arrays.asList(SampleNameListString.split(","));
		String[] output_fields = {"chr", "start", "end", "interval_id", "num_peaks", "num_samples"};
		Stream<String> whole_stream = Arrays.stream(output_fields);
		whole_stream  = Stream.concat(whole_stream, Arrays.stream(SampleNameListString.split(",")).map( n -> n + ".bool" ));
		whole_stream  = Stream.concat(whole_stream, Arrays.stream(SampleNameListString.split(",")).map( n -> n + ".fc"	 ));
		whole_stream  = Stream.concat(whole_stream, Arrays.stream(SampleNameListString.split(",")).map( n -> n + ".qval" ));
		whole_stream  = Stream.concat(whole_stream, Arrays.stream(SampleNameListString.split(",")).map( n -> n + ".pval" ));
		whole_stream  = Stream.concat(whole_stream, Arrays.stream(SampleNameListString.split(",")).map( n -> n + ".start"));
		whole_stream  = Stream.concat(whole_stream, Arrays.stream(SampleNameListString.split(",")).map( n -> n + ".end"  ));
		whole_stream  = Stream.concat(whole_stream, Arrays.stream(SampleNameListString.split(",")).map( n -> n + ".summit"  ));
		output_fields = whole_stream.toArray(String[]::new);
			
		try {
			Scanner 	MergedIntervalFileReader 	= new Scanner(MergedIntervalFile);
			FileWriter 	OutFileWriter 				= new FileWriter(OutFile);
			OutFileWriter.write(String.join("\t", output_fields)+"\n");
			
			String[] 	SummitFiles 				= new String[SampleNameList.size()];
			for(int i = 0; i < SampleNameList.size(); i++){SummitFiles[i]=SampleNameList.get(i)+".summit";}
				
			while(MergedIntervalFileReader.hasNextLine()) {
					
				// READ COLUMN VALUES
				String MergedIntervalString = MergedIntervalFileReader.nextLine();
				String[] MergedInterval = MergedIntervalString.strip().split("\t");
				String 	 chromID		= MergedInterval[0];
				int 	 mstart 		= Integer.parseInt(MergedInterval[1]);
				int 	 mend 			= Integer.parseInt(MergedInterval[2]);
				int[] 	 starts 		= Arrays .stream  (MergedInterval[3] .split(",")).mapToInt		(Integer::parseInt).toArray();
				int[] 	 ends 			= Arrays .stream  (MergedInterval[4] .split(",")).mapToInt		(Integer::parseInt).toArray();
				String[] names 			= 				   MergedInterval[5] .split(",");
				double[] fcs 			= Arrays .stream  (MergedInterval[8] .split(",")).mapToDouble	(Double::parseDouble).toArray() ;
				double[] pvals 			= Arrays .stream  (MergedInterval[9] .split(",")).mapToDouble	(Double::parseDouble).toArray() ;
				double[] qvals 			= Arrays .stream  (MergedInterval[10].split(",")).mapToDouble	(Double::parseDouble).toArray() ;
				int[] 	 summits		= (isNarrow)?
												Arrays .stream(MergedInterval[11].split(",")).mapToInt(Integer::parseInt).toArray():null;
					
				// GROUPS SAMPLES BY EXPERIMENTAL CONDITION
				HashMap<String, ArrayList<String>> groupDict = new HashMap<String, ArrayList<String>>();
				String[] sample_ids = Arrays.stream(names).map(n -> 
																{
																String[] tokens = n.split("_");
																String[] new_tokens = IntStream.range(0,1).
																						mapToObj(i -> tokens[i]).
																							toArray(String[]::new);
																String new_name = String.join("_", new_tokens);
																return new_name;
																}).toArray(String[]::new);
					
				for(String sample_id: sample_ids) {
					String[] tokens_sample_id = sample_id.split("_"); 
					String[] tokens_group_id = IntStream.range(0,1).mapToObj(i -> tokens_sample_id[i]).toArray(String[]::new);
					String group_id = String.join("_", tokens_group_id);
						
					if(!groupDict.containsKey(group_id)) {
						ArrayList<String> value = new ArrayList<String>();
						value.add(sample_id);
						groupDict.put(group_id, value);
					} else {
						groupDict.get(group_id).add(sample_id);
					}
				}
					
				// GET SAMPLES THAT PASS REPLICATE THRESHOLD
				ArrayList<String> passRepThreshList = new ArrayList<String>();
				for(String key_of_replicates: groupDict.keySet()) {
					if(groupDict.get(key_of_replicates).size() >= minRep) {
						passRepThreshList.addAll(groupDict.get(key_of_replicates));
					}
				}
					
				// GET VALUES FROM INDIVIDUAL PEAK SETS
				HashMap<String, ArrayList<String>> 	fcDict 		= new HashMap<String, ArrayList<String>>(), 
													qvalDict 	= new HashMap<String, ArrayList<String>>(), 
													pvalDict 	= new HashMap<String, ArrayList<String>>(), 
													startDict 	= new HashMap<String, ArrayList<String>>(), 
													endDict 	= new HashMap<String, ArrayList<String>>(), 
													summitDict 	= new HashMap<String, ArrayList<String>>();
					
				IntStream.range(0, names.length).forEach(idx -> {
						if(passRepThreshList.contains(sample_ids[idx])) {
							if(fcDict.get(sample_ids[idx]) == null){
								fcDict.put(sample_ids[idx], new ArrayList<String>());
							}
							fcDict.get(sample_ids[idx]).add(String.valueOf(fcs[idx]));
							if(qvalDict.get(sample_ids[idx]) == null){
								qvalDict.put(sample_ids[idx], new ArrayList<String>());
							}
							qvalDict.get(sample_ids[idx]).add(String.valueOf(qvals[idx]));
							if(pvalDict.get(sample_ids[idx]) == null){
								pvalDict.put(sample_ids[idx], new ArrayList<String>());
							}
							pvalDict.get(sample_ids[idx]).add(String.valueOf(pvals[idx]));
							if(startDict.get(sample_ids[idx]) == null){
								startDict.put(sample_ids[idx], new ArrayList<String>());
							}
							startDict.get(sample_ids[idx]).add(String.valueOf(starts[idx]));
							if(endDict.get(sample_ids[idx]) == null){
								endDict.put(sample_ids[idx], new ArrayList<String>());
							}
							endDict.get(sample_ids[idx]).add(String.valueOf(ends[idx]));
							if(summitDict.get(sample_ids[idx]) == null){
								summitDict.put(sample_ids[idx], new ArrayList<String>());
							}
							summitDict.get(sample_ids[idx]).add(String.valueOf(summits[idx]));
						}
					});
						
				// EXTEND COLUMNS SO THAT EACH SAMPLE HAVE A VALUE (NA if none).
				String[] samples = fcDict.keySet().toArray(new String[fcDict.size()]);
				Arrays.sort(samples);
				ArrayList<String> sorted_samples = new ArrayList<String>(Arrays.asList(samples));
					
				ArrayList<String> 	boolList 	= new ArrayList<String>(), 
									fcList 		= new ArrayList<String>(),
									qvalList 	= new ArrayList<String>(), 
									pvalList 	= new ArrayList<String>(), 
									startList 	= new ArrayList<String>(), 
									endList 	= new ArrayList<String>(), 										
									summitList 	= new ArrayList<String>(),
									oList		= new ArrayList<String>();
				
				if(samples.length!=0) {
					int numSamples = samples.length;
					SampleNameList.forEach(
									sample -> {
										if(sorted_samples.contains(sample)) { 
											boolList .add("TRUE"); 
											  fcList .add(String.join(";",    fcDict	.get(sample)));
											qvalList .add(String.join(";",  qvalDict	.get(sample)));
											pvalList .add(String.join(";",  pvalDict	.get(sample)));
										   startList .add(String.join(";", startDict	.get(sample)));
										     endList .add(String.join(";",   endDict	.get(sample)));										   
										  summitList .add(String.join(";",summitDict	.get(sample)));
									} else {
										  boolList 	.add("FALSE")	;	  
										    fcList 	.add("NA")	;	
										  qvalList 	.add("NA")	;	
										  pvalList  .add("NA")	;
										 startList	.add("NA")	;
										   endList  .add("NA")	;	
										summitList 	.add("NA")	;	
									}
								});
						
					oList.add(chromID)		;
					oList.add(mstart+"")	;
					oList.add(mend+"")		;
					oList.add("Interval_"+totalOutIntervals);
					oList.add(names.length+"");
					oList.add(numSamples+"");
					oList.addAll(boolList)	;
					oList.addAll(fcList)	;
					oList.addAll(qvalList)	;
					oList.addAll(pvalList)	;
					oList.addAll(startList)	;
					oList.addAll(endList)	;
					oList.addAll(summitList);
						
					OutFileWriter.write(String.join("\t", oList) + "\n");
					totalOutIntervals++;
				}
				
			}
			MergedIntervalFileReader.close();
			OutFileWriter.close();
		} catch (FileNotFoundException e1) {
				e1.printStackTrace();
		} catch (IOException e) {
				e.printStackTrace();
		}
	}
}

