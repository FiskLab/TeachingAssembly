import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 * 
 * @author Nick 
 * Frag assemble stage 1: It creates the graph structure neccesary
 *         to complete stage 2. It prints out a text representation of the
 *         graph. The graph is constructed based on the minimum overlap
 *         provided. That is, any sequence with a minimum overlap (right of top,
 *         left of bottom) will be added to the other sequences list of
 *         overlapping nodes.
 *         
 *Frag assembly part 2 traverses the graph using a greedy beam search
 *and prints the resulant contigs. 
 * 
 */

// /write a fucntion to see if there are any other things to add
public class FragAssemble {
	public static ArrayList<Contig> Contigs = new ArrayList<Contig>();// will hold contigs
	//public static HashMap<String, Boolean> visited = new HashMap<String, Boolean>();
	public static int minOverlap = 0;
	public static ArrayList<Frag> graph = new ArrayList<Frag>();// the graph (a
																// list of all
																// the parent
																// nodes)
	public static HashMap<String, String> namesNsequences = new HashMap<String, String>();// name
																							// is
																							// suffiecent
																							// to
																							// get
																							// a
																							// sequence--more
	// lightweight than a separate object.
	public static String usage = "Usage: FragAssemble [fragments.fasta][minimum overlap]";

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		if (args.length != 2) {// wrong!
			System.out.println(usage);
		}
		String file = args[0];
		String minOverLapString = args[1];
		try {
			minOverlap = Integer.parseInt(minOverLapString);
		} catch (NumberFormatException e) {
			System.err.println("invalid minimum overlap");
			System.exit(1);
		}
		minOverlap-=1;
		BufferedReader read = new BufferedReader(new FileReader(file));
		String line = null;
		String[] lineArray = null;
		String tempSeq = "";
		String nameTemp = "";
		// read the file in line by line (well, more like 2 lines at a time)
		while ((line = read.readLine()) != null) {
			if (line.length() > 2) {
				if (line.charAt(0) == '>') {
					line = line.trim();
					nameTemp = line.substring(1);
					line = read.readLine();
					tempSeq = line.trim();
					if (checkIfValidSeq(tempSeq) == true) {// make sure these
															// are good
															// nucleotide
															// sequences
						namesNsequences.put(nameTemp, tempSeq);// stick them in
																// the hashmap
						// printmap.put(nameTemp, false); an artifact from an
						// old printing problem.
					}
				}
			}
		}
		Frag tempFrag = null;
		for (String parent : namesNsequences.keySet()) {// every sequence will
														// become a parent!
			tempFrag = new Frag(parent, namesNsequences.get(parent));// this
																		// will
																		// be
																		// that
																		// parent
			String parSeq = "";
			String chiSeq = "";
			for (String potentialChild : namesNsequences.keySet()) {// gotta
																	// check all
																	// the other
																	// sequences
				if (!((potentialChild).equals(parent))) {// otherwise, itll
															// always overlap
															// with itself XD
					// store the sequences
					parSeq = namesNsequences.get(parent);
					chiSeq = namesNsequences.get(potentialChild);
					// checking to see if there is any overlap. I may revise
					// this
					// to start at the end of the parent sequence and work
					// backwards (this will
					// really cut down on computation time).
					String check = "";
					String bestOverlap = "";
					String potentialOverlap = "";
					String potentialOverlapChi = "";
					for (int i = 0; i < parSeq.length(); i++) {
						potentialOverlap = parSeq.substring(i);
						potentialOverlapChi = "";
						if (potentialOverlap.length() <= chiSeq.length()) {
							potentialOverlapChi = chiSeq.substring(0,
									potentialOverlap.length());
						} else {
							potentialOverlapChi = chiSeq;
						}
						check = parSeq
								.substring((parSeq.length() - potentialOverlap
										.length()));
						if (check.equals(potentialOverlap)
								&& potentialOverlapChi.equals(potentialOverlap)) {
							if (potentialOverlap.length() > bestOverlap
									.length()) {
								bestOverlap = potentialOverlap;
							}
						}
					}
					// are you good enough for our user?
					// if so, the child can be added to the parent.
					if (bestOverlap.length() > minOverlap) {
						Frag child = new Frag(potentialChild,
								namesNsequences.get(potentialChild));
						child.overlapSeq = bestOverlap;
						child.amountOverlap = bestOverlap.length();
						tempFrag.overlapFrags.add(child);
					}
				}
			}
			// add the parent to the overall structure
			graph.add(tempFrag);

		}
		//////////////////////START PART II///////////////////////////////
		int conSize = 0;
		int conSizeNew = 0;
		while (true) {//have a break statement when the number
			///of contigs remain constant down below.
			while (anythingLeft()) {//while at least one node has a child
				beamSearch();//perform the beamSearch
			}
			addUnusedFrags(); //add unused fragments to contig list
			conSizeNew = Contigs.size(); 
			if (conSizeNew == conSize) {
				break; //if it is no fewer than orignal, break
			} else {
				conSize = conSizeNew;
			}
			
			rebuildGraph(); //rebuilds the graph from contigs
			Contigs.clear();//clears the current contig list
		}
		
		addUnusedFrags(); //add the last round of unused fragments to contigs
		int contigCount = 1; //for printing
		///print all the contigs
		for (Contig c : Contigs) {
			System.out.println("Contig: " + contigCount);
			System.out.println("\t Consensus:");
			int lineStart = 0;
			int lineEnd = 0;
			while (lineEnd < c.consensus.length()) {
				lineEnd += 45;
				if (lineEnd > c.consensus.length()) {
					lineEnd = c.consensus.length();
				}
				System.out.println("\t \t"
						+ c.consensus.substring(lineStart, lineEnd));
				lineStart = lineEnd;
			}
			// System.out.println("\t Consensus: "+c.consensus);
			System.out.println("\t Size of Contig: " + c.consensus.length());
			
			if(c.frags.size()<=3){
				System.out.println("\t Made of Fragments " + c.frags+"\n");
			}
			else{
				System.out.println("\t Made of Fragments:");
				System.out.print("\t \t [");
				int count=0;
				for(String f: c.frags){
					System.out.print(f);
					if(count>8){
						System.out.println();
						System.out.print("\t \t");
						count=0;
					}
				}
				System.out.println(" ]");
				
			}
			contigCount++;
			// }
		}
	}

	/*
	 * Checks if a given sequence is a valid DNA sequence.
	 */
	public static boolean checkIfValidSeq(String seq) {
		boolean holdMyBoolSon = true;
		for (char c : seq.toCharArray()) {
			if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
				System.out
						.println("Invalid nucleotide in sequence. Sequence will be ignored");
				holdMyBoolSon = false;
			}
		}
		return holdMyBoolSon;

	}
	/*
	 * 
	 * Checks to see if there is at least one parent
	 * left that has a child.
	 */
	public static boolean anythingLeft() {
		for (Frag f : graph) {
			if (!f.overlapFrags.isEmpty()) {
				return true;
			}
		}

		return false;
	}

	/*
	 * Prints the graph, used in part one. 
	 */
	public static void printGraph() {
		String toPrint = "";
		for (Frag seq : graph) {
			toPrint = "SeqID:" + seq.name + " has overlaps with [";
			if (seq.overlapFrags.size() == 0) {
				toPrint += "<none>]";
			} else {
				for (Frag child : seq.overlapFrags) {
					if (child != seq) {
						toPrint += ("(" + "SeqName" + child.name
								+ "| OverlapSeq " + child.overlapSeq + " ),");
					}
				}
				toPrint += "]";
			}
			System.out.println(toPrint);
			toPrint = "";
		}
	}

	/*
	 * condenses a child and a parent into just the parent
	 * The name of the resultant Frag will be that of the parent
	 * and the children of the resultant frag will be that of the 
	 * children. This way, anything that has the parent as a child
	 * will still recognize this condensed fragment, but the children
	 * are appropriate. 
	 */
	public static Frag condenseFrags(Frag top, Frag bottom, String overlapSeq) {
		String topSeq = top.seq;
		String bottomSeq = bottom.seq;
		int overlapLen = overlapSeq.length();
		topSeq = topSeq.substring(0, topSeq.length() - overlapLen);
		top.seq = topSeq + bottomSeq;
		top.overlapFrags = bottom.overlapFrags;
		return top;
	}

	/*
	 * Simply finds the fragment that has the highest overlap or edge score. 
	 */
	public static Frag findMaxFrag() {
		int max = 0;
		Frag first = new Frag("", "");
		for (Frag f : graph) {
			if (f.overlapFrags.size() > 0) {
				for (Frag child : f.overlapFrags) {
					if (child.amountOverlap > max) {
						max = child.amountOverlap;
						first = f;
					}
				}
			}
		}

		return first;
	}

	/*
	 * findBestOverlap finds the child with the maximum overlap with the given parent
	 * from its defined list of frags. In other words, finds the best edge
	 * for a given node. 
	 */
	public static Frag findBestOverlap(Frag frag) {
		Frag bestOver = new Frag("", "");
		Frag tempOver = new Frag("", "");
		int curMax = 0;
		for (Frag f : frag.overlapFrags) {
			if (f.amountOverlap > curMax) {
				curMax = f.amountOverlap;
				tempOver = f;
				for (int i = 0; i < graph.size(); i++) {
					if (graph.get(i).name.equals(tempOver.name)) {
						bestOver = graph.get(i);
					}
				}
			}
		}

		return bestOver;
	}

	/*
	 * beamSearch performs the main functionality of the assembly.
	 * It finds the best fragment and the fragment that best overlaps with
	 * it, joins them together, removes the child from the list, sets the childs
	 * children to the joint node (which keeps the parent's name) and continues
	 * joining until there are no more children in this chain. 
	 */
	public static void beamSearch() {
		Frag best = new Frag("", "");
		Frag overlap = new Frag("", "");
		best = findMaxFrag();//node with highest edge
		Contig contig = new Contig();//new contig to hold this info
		contig.frags.add(best.name);//the best automatically gets added
		overlap = findBestOverlap(best);//find its best overlap
		//////this is the case if there are no children//////
		if (overlap.name.equals("")) {
			contig.consensus = best.seq;
			for (Iterator<Frag> iterator2 = graph.iterator(); iterator2
					.hasNext();) {
				Frag f = iterator2.next();
				if (f.name.equals(best.name)) {
					iterator2.remove();
				}
			}
		}
		////////This is the case if there is at least one child////////
		while (!overlap.name.equals("") && !overlap.seq.equals("")) {//until no children left
			contig.frags.add(overlap.name);
			Frag newFrag = condenseFrags(best, overlap, overlap.finalOverlap);
			for (Iterator<Frag> iterator = graph.iterator(); iterator.hasNext();) {
				Frag f = iterator.next();
				if (f.name.equals(overlap.name)) {
					iterator.remove();
				}
				if (f.name.equals(best.name)) {
					iterator.remove();
				}
			}
			best = newFrag;
			overlap = findBestOverlap(best);
		}
		contig.consensus = best.seq;
		Contigs.add(contig);
		
	}

	/**
	 * Take the frags that weren't matched with anything in a given
	 * beamSearch iteration and add them to our list of contigs (
	 * as the lab requires this)
	 */
	public static void addUnusedFrags() {
		Contig contig = new Contig();
		for (Frag f : graph) {
			contig = new Contig();
			contig.consensus = f.seq;
			contig.frags.add(f.name);
			Contigs.add(contig);
		}
		graph.clear();

	}
	/*
	 * Rebuild graph does the same thing that main does when 
	 * reading in the file, except it uses all the sequences 
	 * in the list of contigs to do so. It reads in the contig
	 * names and sequences, makes Frag s out of them, and makes 
	 * makes a graph out of them.
	 */
	public static void rebuildGraph() {
		namesNsequences.clear();
		for (Contig c : Contigs) {
			namesNsequences.put(
					c.frags.toString().substring(1,
							c.frags.toString().length() - 1), c.consensus);
		}
		Frag tempFrag = null;
		for (String parent : namesNsequences.keySet()) {// every sequence will
														// become a parent!
			tempFrag = new Frag(parent, namesNsequences.get(parent));// this
																		// will
																		// be
																		// that
																		// parent
			String parSeq = "";
			String chiSeq = "";
			for (String potentialChild : namesNsequences.keySet()) {// gotta
																	// check all
																	// the other
																	// sequences
				if (!((potentialChild).equals(parent))) {// otherwise, itll
															// always overlap
															// with itself XD
					// store the sequences
					parSeq = namesNsequences.get(parent);
					chiSeq = namesNsequences.get(potentialChild);
					// checking to see if there is any overlap. I may revise
					// this
					// to start at the end of the parent sequence and work
					// backwards (this will
					// really cut down on computation time).
					String check = "";
					String bestOverlap = "";
					String potentialOverlap = "";
					String potentialOverlapChi = "";
					for (int i = 0; i < parSeq.length(); i++) {
						potentialOverlap = parSeq.substring(i);
						potentialOverlapChi = "";
						if (potentialOverlap.length() <= chiSeq.length()) {
							potentialOverlapChi = chiSeq.substring(0,
									potentialOverlap.length());
						} else {
							potentialOverlapChi = chiSeq;
						}
						check = parSeq
								.substring((parSeq.length() - potentialOverlap
										.length()));
						if (check.equals(potentialOverlap)
								&& potentialOverlapChi.equals(potentialOverlap)) {
							if (potentialOverlap.length() > bestOverlap
									.length()) {
								bestOverlap = potentialOverlap;
							}
						}
					}
					// are you good enough for our user?
					// if so, the child can be added to the parent.
					if (bestOverlap.length() > minOverlap) {
						Frag child = new Frag(potentialChild,
								namesNsequences.get(potentialChild));
						child.overlapSeq = bestOverlap;
						child.amountOverlap = bestOverlap.length();
						tempFrag.overlapFrags.add(child);
					}
				}
			}
			// add the parent to the overall structure
			graph.add(tempFrag);

		}

	}
}
