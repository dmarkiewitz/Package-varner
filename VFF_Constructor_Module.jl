module VFF_Constructor_Module
export metabolite
export reaction_function
export VFF_constructor_function
export add_reaction
export delete_reaction
export adjust_reaction
export Reaction_Gene_Map_function
export Gene_Symbol_Order_function
export Rule_building
export Rule_edit
export VFF_to_mat

using DataFrames
using CSV
using LibCURL
using DelimitedFiles
using MAT

	#builds the information needed for the metabolite section of a VFF File 
  @doc """
	metabolite(File_Name_txt::String,Decode_File_txt::String,rxns_column_number::Int,first_rxn_row::Int,last_rxn_row::Int,rxn_ref_name_column::Int,lowerbound_reaction_column::Int,upperbound_reaction_column::Int,Sym_colm::Int,Name_colm::Int) \n
	File_Name_txt-> Where the reactions are saved in a tab deliminated .txt file \n
	Decode_File_txt-> Where the file where the chemical compounds and their symbolic representation is list, this should be a tab deliminated .txt file \n
	rxns_column_number-> The column in which the symbolic representation of the reaction is located \n
	first_rxn_row-> The first row containing a reaction in the .txt file \n 
	last_rxn_row-> The last row containing a reaction in the .txt file \n
	rxn_ref_name_column-> The column number which contains the reaction name column \n
	lowerbound_reaction_column-> The column number contianing the lowerbound value for the flux through an arbitrary reaction \n
	upperbound_reaction_column-> The column number contianing the upperbound value for the flux through an arbitrary reaction \n
	Sym_colm-> The column number contianing the symbolic name \n
	Name_colm-> The column number which contians the full name used for decoding the symbolic name \n

	This function produces a string array which has the metabolite information in .VFF format.

       """ ->
	function metabolite(File_Name_txt::String,Decode_File_txt::String,rxns_column_number::Int,first_rxn_row::Int,last_rxn_row::Int,rxn_ref_name_column::Int,
	lowerbound_reaction_column::Int,upperbound_reaction_column::Int,Sym_colm::Int,Name_colm::Int)
	    
	    #Collection of data
	    rxn_path_data=readdlm(File_Name_txt, '\t', String, '\n')

	    #seperation of data
	    rxn_names=rxn_path_data[first_rxn_row:last_rxn_row,rxn_ref_name_column]#finds how rxn are named used for reference
	    rxns=rxn_path_data[first_rxn_row:last_rxn_row,rxns_column_number]#the actual stoichometric reactions

	    #building unique metabolic list
	    chemicals_raw=[]
	    for i in 1:length(rxns)#length(rxns) #column i of matrix
		bigsplit=split(rxns[i],r"  -> |  <=> ")#splits rxn into left and right half
		left_rxn=bigsplit[1]
		left_chem_raw=split(left_rxn," + ")#subdivides the left half into chemicals and their coeffiecent
		left_chem_and_co=split.(left_chem_raw," ")#splits up chemicals and their coeffiecent
		for j in 1:length(left_chem_and_co)#collects all unique chemicals on left side
		    if length(left_chem_and_co[j])==2 && (!(in(left_chem_and_co[j][2],chemicals_raw)) && left_chem_and_co[j][2]!="")#if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
		        chemicals_raw=push!(chemicals_raw,left_chem_and_co[j][2])
		    elseif length(left_chem_and_co[j])==1 && (!(in(left_chem_and_co[j][1],chemicals_raw)) && left_chem_and_co[j][1]!="")#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
		        chemicals_raw=push!(chemicals_raw,left_chem_and_co[j][1])
		    end
		end
		right_rxn=bigsplit[2]
		right_chem_raw=split(right_rxn," + ")#subdivides the right half into chemicals and their coeffiecent
		right_chem_and_co=split.(right_chem_raw)#splits up chemicals and their coeffiecent
		for j in 1:length(right_chem_and_co)#collects all unique chemicals on right side
		    if length(right_chem_and_co[j])==2 && !(in(right_chem_and_co[j][2],chemicals_raw)) && !(right_chem_and_co[j][2]=="")#if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
		        chemicals_raw=push!(chemicals_raw,right_chem_and_co[j][2])
		    elseif length(right_chem_and_co[j])==1 && !(in(right_chem_and_co[j][1],chemicals_raw)) && !(right_chem_and_co[j][1]=="")#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
		        chemicals_raw=push!(chemicals_raw,right_chem_and_co[j][1])
		    end
		end
	    end

	    #prebuilding of chemical's compartment organization vectors
	    c_unsort=[]
	    m_unsort=[]
	    r_unsort=[]
	    x_unsort=[]
	    l_unsort=[]
	    e_unsort=[]

	    for i in 1:length(chemicals_raw)#length(rxns) #column i of matrix
		#Bin sorting of chemicals by their compartmental location
		if occursin("[c]",chemicals_raw[i])
		    c_unsort=push!(c_unsort,chemicals_raw[i])
		elseif occursin("[m]",chemicals_raw[i])
		    m_unsort=push!(m_unsort,chemicals_raw[i])
		elseif occursin("[r]",chemicals_raw[i])
		    r_unsort=push!(r_unsort,chemicals_raw[i])
		elseif occursin("[x]",chemicals_raw[i])
		    x_unsort=push!(x_unsort,chemicals_raw[i])
		elseif occursin("[l]",chemicals_raw[i])
		    l_unsort=push!(l_unsort,chemicals_raw[i])
		elseif occursin("[e]",chemicals_raw[i])
		    e_unsort=push!(e_unsort,chemicals_raw[i])
		end
	    end

	    #Sorting chemicals via ASCII assignments
	    c_sorted=sort(c_unsort)
	    m_sorted=sort(m_unsort)
	    r_sorted=sort(r_unsort)
	    x_sorted=sort(x_unsort)
	    l_sorted=sort(l_unsort)
	    e_sorted=sort(e_unsort)

	    #Prebuilding sorted chemical vector
	    chemicals_sorted=[]

	    #Building sorted chemical vector by adding blocks of sorted chemicals
	    chemicals_sorted=c_sorted
	    chemicals_sorted=append!(chemicals_sorted,m_sorted)
	    chemicals_sorted=append!(chemicals_sorted,r_sorted)
	    chemicals_sorted=append!(chemicals_sorted,x_sorted)
	    chemicals_sorted=append!(chemicals_sorted,l_sorted)
	    chemicals_sorted=append!(chemicals_sorted,e_sorted)

	    #creation of V_boundaries_vector
	    V_boundaries_vector=parse.(Float64,rxn_path_data[first_rxn_row:last_rxn_row,lowerbound_reaction_column])#column one is lower bound and column two is upper bound,rows correspond to rxn_names ordering
	    V_boundaries_vector=cat(V_boundaries_vector,(parse.(Float64,rxn_path_data[first_rxn_row:last_rxn_row,upperbound_reaction_column])); dims=2)

	    #collecting the decoding document
	    decoded_data=readdlm(Decode_File_txt, '\t', String, '\n')
	    ddrs,ddcs=size(decoded_data)
	    
	    #Seperating out data
	    coded_data=decoded_data[2:ddrs,Sym_colm]
	    decoding_info=decoded_data[2:ddrs,Name_colm]

	    #Prebuilding (symbolic name, real name, kegg ID)
	    Pre_VFF_metabolic=Array{Union{Nothing, String}}(nothing, length(chemicals_sorted), 3)
	    

	    #Note bug not yet fixed was certain chemicals in the decoding are not the proper or identifiable name e.g. N-(omega)-Hydroxyarginine in the decoding and its actual name N(omega)-Hydroxyarginine
	    for g in 1:length(chemicals_sorted)
		#sorting of chemicals by their compartmental location
		#decoding the symboles and then quering kegg
		Pre_VFF_metabolic[g,1]=chemicals_sorted[g]
		Pre_VFF_metabolic[g,2]="nothing"
		Pre_VFF_metabolic[g,3]="nothing"
		for h in 1:length(coded_data)
		    if isequal(chemicals_sorted[g],coded_data[h]) ### correction from previous version
		        if occursin("[c]",chemicals_sorted[g])
		            real_name=split(decoding_info[h]," cytosol")[1]
		            Pre_VFF_metabolic[g,2]=real_name
			    real_name=replace(real_name, " " => "%20") ###correction from previous version
		            h=read(`curl -g http://rest.kegg.jp/find/compound/$real_name/`, String)
		            name=split(h,"\t")
				if length(name)>=2
				    Cid=name[1]
				    Pre_VFF_metabolic[g,3]=Cid
				end
		        elseif occursin("[m]",chemicals_sorted[g])
		            real_name=split(decoding_info[h]," mitochondria")[1]
		            Pre_VFF_metabolic[g,2]=real_name
			    real_name=replace(real_name, " " => "%20") ###correction from previous version
		            h=read(`curl -g http://rest.kegg.jp/find/compound/$real_name/`, String)
		            name=split(h,"\t")
				if length(name)>=2
				    Cid=name[1]
				    Pre_VFF_metabolic[g,3]=Cid
				end
		        elseif occursin("[r]",chemicals_sorted[g])
		            real_name=split(decoding_info[h]," ER")[1]
		            Pre_VFF_metabolic[g,2]=real_name
			    real_name=replace(real_name, " " => "%20") ###correction from previous version
		            h=read(`curl -g http://rest.kegg.jp/find/compound/$real_name/`, String)
		            name=split(h,"\t")
				if length(name)>=2
				    Cid=name[1]
				    Pre_VFF_metabolic[g,3]=Cid
				end
		        elseif occursin("[x]",chemicals_sorted[g])
		            real_name=split(decoding_info[h]," peroxisome")[1]
		            Pre_VFF_metabolic[g,2]=real_name
			    real_name=replace(real_name, " " => "%20") ###correction from previous version
		            h=read(`curl -g http://rest.kegg.jp/find/compound/$real_name/`, String)
		            name=split(h,"\t")
				if length(name)>=2
				    Cid=name[1]
				    Pre_VFF_metabolic[g,3]=Cid
				end
		        elseif occursin("[l]",chemicals_sorted[g])
		            real_name=split(decoding_info[h]," lysome")[1]
		            Pre_VFF_metabolic[g,2]=real_name
			    real_name=replace(real_name, " " => "%20") ###correction from previous version
		            h=read(`curl -g http://rest.kegg.jp/find/compound/$real_name/`, String)
		            name=split(h,"\t")
				if length(name)>=2
				    Cid=name[1]
				    Pre_VFF_metabolic[g,3]=Cid
				end
		        elseif occursin("[e]",chemicals_sorted[g])
		            real_name=split(decoding_info[h]," extracellular")[1]
		            Pre_VFF_metabolic[g,2]=real_name
			    real_name=replace(real_name, " " => "%20") ###correction from previous version
		            h=read(`curl -g http://rest.kegg.jp/find/compound/$real_name/`, String)
		            name=split(h,"\t")
				if length(name)>=2
				    Cid=name[1]
				    Pre_VFF_metabolic[g,3]=Cid
				end
		        end
		    end
		end
	    end

		#Pre-construting VFF array to be saved to external file which can be later read and then constructed into the complete .VFF for the corresponding model
		VFF_metabolite_array=Array{Union{Nothing, String}}(nothing, length(Pre_VFF_metabolic[:,1]), 1)

		for b in 1:length(Pre_VFF_metabolic[:,1])

			#assigning symbolic name
			sym_name=Pre_VFF_metabolic[b,1]
			
			#assigning real name
			real_name=Pre_VFF_metabolic[b,2]	

			#assigning Kegg ID
			Kegg_ID=Pre_VFF_metabolic[b,3]

			#Concatenating into correct .VFF format
			VFF_metabolite="$sym_name"*"="*"$real_name"*"::"*"$Kegg_ID"

			#Adding it to the VFF_metabolite_array
			VFF_metabolite_array[b]=VFF_metabolite
		end

		return VFF_metabolite_array
	end


	#builds the information needed for the reaction section of a VFF File 
  @doc """
	reaction_function(File_Name_txt::String,rxns_column_number::Int,first_rxn_row::Int,last_rxn_row::Int,rxn_ref_name_column::Int,lowerbound_reaction_column::Int,upperbound_reaction_column::Int,Descriptor_column::Int) \n
	File_Name_txt-> Where the reactions are saved in a tab deliminated .txt file \n
	rxns_column_number-> The column in which the symbolic representation of the reaction is located \n
	first_rxn_row-> The first row containing a reaction in the .txt file \n
	last_rxn_row-> The last row containing a reaction in the .txt file \n
	rxn_ref_name_column-> The column number which contains the reaction name column \n
	lowerbound_reaction_column-> The column number contianing the lowerbound value for the flux through an arbitrary reaction \n
	upperbound_reaction_column-> The column number contianing the upperbound value for the flux through an arbitrary reaction \n
	Descriptor_column-> The column number which contains the reaction description which can be used to query for e.c's \n

	This function produces a string array with the information about the reactions in the network of reaction in .VFF format.

       """ ->
	function reaction_function(File_Name_txt::String,rxns_column_number::Int,first_rxn_row::Int,last_rxn_row::Int,rxn_ref_name_column::Int,
		                    lowerbound_reaction_column::Int,upperbound_reaction_column::Int,Descriptor_column::Int)

	    #Collection of data
	    rxn_path_data=readdlm(File_Name_txt,'\t',String,'\n')

	    #seperation of data
	    rxn_names=rxn_path_data[first_rxn_row:last_rxn_row,rxn_ref_name_column]#finds how rxn are named used for reference
	    rxns=rxn_path_data[first_rxn_row:last_rxn_row,rxns_column_number]#the actual stoichometric reactions

	    #building reactions list
	    reactants_list=Array{Union{Nothing, String}}(nothing, length(rxns), 1)
	    products_list=Array{Union{Nothing, String}}(nothing, length(rxns), 1)
	    for i in 1:length(rxns)#length(rxns) #column i of matrix
		bigsplit=split(rxns[i],r"  -> |  <=> ")#splits rxn into left and right half
		left_rxn=bigsplit[1]
		left_chem_raw=split(left_rxn," + ")#subdivides the left half into chemicals and their coeffiecent
		left_chem_and_co=split.(left_chem_raw," ")#splits up chemicals and their coeffiecent
		reactants=""
		products=""
		for j in 1:length(left_chem_and_co)#collects all unique chemicals on left side
			if j!=length(left_chem_and_co)       
			    if length(left_chem_and_co[j])==2 && left_chem_and_co[j][2]!="" #if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
				coeff=left_chem_and_co[j][1]
				chemical=left_chem_and_co[j][2]			
				reactants=reactants*"$coeff"*"*"*"$chemical"*"+"
			    elseif length(left_chem_and_co[j])==1 && left_chem_and_co[j][1]!=""#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
				chemical=left_chem_and_co[j][1]			
				reactants=reactants*"$chemical"*"+"
			    end
			else
			    if length(left_chem_and_co[j])==2 && left_chem_and_co[j][2]!="" #if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
				coeff=left_chem_and_co[j][1]
				chemical=left_chem_and_co[j][2]			
				reactants=reactants*"$coeff"*"*"*"$chemical"
			    elseif length(left_chem_and_co[j])==1 && left_chem_and_co[j][1]!=""#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
				chemical=left_chem_and_co[j][1]			
				reactants=reactants*"$chemical"
			    end
			end
		end
		reactants_list[i]=reactants
		right_rxn=bigsplit[2]
		right_chem_raw=split(right_rxn," + ")#subdivides the right half into chemicals and their coeffiecent
		right_chem_and_co=split.(right_chem_raw)#splits up chemicals and their coeffiecent
		for j in 1:length(right_chem_and_co)#collects all unique chemicals on right side
			if j!=length(right_chem_and_co)          
			   if length(right_chem_and_co[j])==2 && right_chem_and_co[j][2]!=""#if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
				coeff=right_chem_and_co[j][1]
				chemical=right_chem_and_co[j][2]			
				products=products*"$coeff"*"*"*"$chemical+"
			    elseif length(right_chem_and_co[j])==1 && right_chem_and_co[j][1]!=""#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
				chemical=right_chem_and_co[j][1]			
				products=products*"$chemical+"
			    end
			else
			    if length(right_chem_and_co[j])==2 && right_chem_and_co[j][2]!=""#if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
				coeff=right_chem_and_co[j][1]
				chemical=right_chem_and_co[j][2]			
				products=products*"$coeff"*"*"*"$chemical"
			    elseif length(right_chem_and_co[j])==1 && right_chem_and_co[j][1]!=""#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
				chemical=right_chem_and_co[j][1]			
				products=products*"$chemical"
			    end
			end	
		end
		products_list[i]=products
	    end

	    #creation of V_boundaries_vector
	    V_boundaries_vector=parse.(Float64,rxn_path_data[first_rxn_row:last_rxn_row,lowerbound_reaction_column])#column one is lower bound and column two is upper bound,rows correspond to rxn_names ordering
	    V_boundaries_vector=cat(V_boundaries_vector,(parse.(Float64,rxn_path_data[first_rxn_row:last_rxn_row,upperbound_reaction_column])); dims=2)

	    #Seperating up data
	    Descriptor_list=rxn_path_data[first_rxn_row:last_rxn_row,Descriptor_column]#finds how rxns enzymes are described for later quering

	    #Initializing E.C. Storage List
	    ec_list=Array{Union{Nothing, String}}(nothing, length(Descriptor_list), 1)

	    #obtaining E.C list
	    for v in 1:length(Descriptor_list)
		ec_temp=[]
		Descriptor=Descriptor_list[v]
	    	pos_ecs=read(`curl -X GET http://rest.kegg.jp/find/enzyme/$Descriptor/`, String)
		sep_pos_ecs=split(pos_ecs,"\n")
			for u in 1:length(sep_pos_ecs)
				s_ec=split(sep_pos_ecs[u],"\t")
			   if length(s_ec)!=1
				r_descriptors=split(s_ec[2],"; ")#retrived discriptors of enzymes
				for w in 1:length(r_descriptors)
					if isequal(Descriptor,r_descriptors[w])
	 					ec_temp=push!(ec_temp,s_ec[1])
					end
				end
			   end
			end
		ec_for_rxn=""	
		for s in 1:length(ec_temp)
			temp=ec_temp[s]		
			ec_for_rxn="$ec_for_rxn"*"$temp"*"::"
		end
		ec_for_rxn=ec_for_rxn*"[]"
		ec_list[v]=ec_for_rxn	
	    end
	    
	    #Initializing reaction .vff Storage List
	    rxn_vff_list=Array{Union{Nothing, String}}(nothing, length(Descriptor_list), 1)

		
	    #constructing reaction .vff string
	    for t in 1:length(rxn_names)
		rxn_temp=rxn_names[t]
		ec_temp=ec_list[t]
		reactant_temp=reactants_list[t]
		product_temp=products_list[t]
		lb_temp=V_boundaries_vector[t,1]
		ub_temp=V_boundaries_vector[t,2]

		#constructing full string
		rxn_vff_temp="$rxn_temp,"*"$ec_temp,"*"$reactant_temp,"*"$product_temp,"*"$lb_temp,"*"$ub_temp"
		
		#saving constructed string
		rxn_vff_list[t]=rxn_vff_temp
	    end
	    return rxn_vff_list #returns the reaction information as a vector of strings in vff format such that saving is the only required action
	end


  @doc """	
	VFF_constructor_function(File_Name_txt::String,Decode_File_txt::String,rxns_column_number::Int,first_rxn_row::Int,last_rxn_row::Int,rxn_ref_name_column::Int,lowerbound_reaction_column::Int,upperbound_reaction_column::Int,Sym_colm::Int,Name_colm::Int,Descriptor_column::Int,KOC::String,Rule_Document::String,Rule_String_Array::Array,Save_VFF_name::String) \n
	File_Name_txt-> Where the reactions are saved in a tab deliminated .txt file \n
	Decode_File_txt-> Where the file where the chemical compounds and their symbolic representation is list, this should be a tab deliminated .txt file \n
	rxns_column_number-> The column in which the symbolic representation of the reaction is located \n
	first_rxn_row-> The first row containing a reaction in the .txt file \n
	last_rxn_row-> The last row containing a reaction in the .txt file \n
	rxn_ref_name_column-> The column number which contains the reaction name column \n
	lowerbound_reaction_column-> The column number contianing the lowerbound value for the flux through an arbitrary reaction \n
	upperbound_reaction_column-> The column number contianing the upperbound value for the flux through an arbitrary reaction \n
	Sym_colm-> The column number contianing the symbolic name \n
	Name_colm-> The column number which contians the full name used for decoding the symbolic name \n
	Descriptor_column-> The column number which contains the reaction description which can be used to query for e.c's \n
	Rule_Document-> The document location / name which contians the information of the rules for the reactions \n
	Rule_String_Array-> A String array that contains the rule information in .VFF format \n
	Save_VFF_name-> The name of the .txt file one wants to create and save to should have .txt in its name \n

	This function constructs the the model in VFF format and then saves it to the Save_VFF_name location where julia is running in.

       """ ->
	function VFF_constructor_function(File_Name_txt::String,Decode_File_txt::String,rxns_column_number::Int,first_rxn_row::Int,last_rxn_row::Int,rxn_ref_name_column::Int,lowerbound_reaction_column::Int,upperbound_reaction_column::Int,Sym_colm::Int,Name_colm::Int,Descriptor_column::Int,KOC::String,Rule_Document::String,Rule_String_Array::Array,Save_VFF_name::String)

	#Constructing VFF format array such that it can be immediately saved as a .txt file
	Temp="// Contents: Chemical reaction string"
	Temp=vcat(Temp,"// Records: Reaction_name,{ecnumbers::[]},reactants,products,reverse,forward")
	Temp=vcat(Temp,"// Records: MDH,ec:1.1.1.37::[],mal_L[c]+nad[c],h[c]+nadh[c]+oaa[c],-inf,inf")
	Temp=vcat(Temp,"// Records: SUCD1m,ec:3.2.1.48::ec:3.2.1.20::ec:3.2.1.3::ec:3.2.1.10::[],q10[m]+succ[m],fum[m]+q10h2[m],0,inf")
	Temp=vcat(Temp,"#REACTION::START --------------------------------------------- //")
	Temp=vcat(Temp,reaction_function(File_Name_txt,rxns_column_number,first_rxn_row,last_rxn_row,rxn_ref_name_column,lowerbound_reaction_column,upperbound_reaction_column,Descriptor_column))
	Temp=vcat(Temp,"#REACTION::END ----------------------------------------------- //")
	Temp=vcat(Temp,"")
	Temp=vcat(Temp,"// Contents: Mapping between metabolite symbol and actual name")
	Temp=vcat(Temp,"// Contents: and KEGG compound ID")
	Temp=vcat(Temp,"// Records: Metabolite_symbol=Metabolite_name::KEGG_ID")
	Temp=vcat(Temp,"// Records: q10[m]=Ubiquinone-10::cpd:C11378")
	Temp=vcat(Temp,"#METABOLITES::START ------------------------------------------ //")

	#Obtaining the required information
	VFF_metabolite_array=metabolite(File_Name_txt,Decode_File_txt,rxns_column_number,first_rxn_row,last_rxn_row,rxn_ref_name_column,lowerbound_reaction_column,upperbound_reaction_column,Sym_colm,Name_colm)
	
	#adding VFF metabolite information to the correct spot
	Temp=vcat(Temp,VFF_metabolite_array)

	Temp=vcat(Temp,"#METABOLITES::END -------------------------------------------- //")
	Temp=vcat(Temp,"")
	Temp=vcat(Temp,"// Contents: Rule records exported from the MAT binary file")
	Temp=vcat(Temp,"// Records: Reaction_name={rule} or [] if none. The rules are")
	Temp=vcat(Temp,"// Records: boolean statements of gene states")
	Temp=vcat(Temp,"// Records: SUCD1m=(x(1481) | x(1918))")
	Temp=vcat(Temp,"// Records: ASPte=[]")
	Temp=vcat(Temp,"#RULES::START ------------------------------------------------ //")
	#Obtaining or creating rules and thus the rule section	
	VFF_Rule_Array=Rule_building(File_Name_txt,rxns_column_number,first_rxn_row,last_rxn_row,rxn_ref_name_column,Rule_Document,Rule_String_Array)
	Temp=vcat(Temp,VFF_Rule_Array)
	Temp=vcat(Temp,"#RULES::END -------------------------------------------------- //")
	Temp=vcat(Temp,"")
	Temp=vcat(Temp,"// Contents: Mapping between reactions in the model and genes")
	Temp=vcat(Temp,"// Records: Reaction_name={:: delimited gene set}::[]")
	Temp=vcat(Temp,"// Records: MTHFC=hsa:4522.1::hsa:80068.1::hsa:285216.1::[]")
	Temp=vcat(Temp,"#REACTION-GENE-MAP::START ------------------------------------ //")
	RGMVFF=Reaction_Gene_Map_function(File_Name_txt,KOC,first_rxn_row,last_rxn_row,rxn_ref_name_column)
	Temp=vcat(Temp,RGMVFF)
	Temp=vcat(Temp,"#REACTION-GENE-MAP::END -------------------------------------- //")
	Temp=vcat(Temp,"")
	Temp=vcat(Temp,"// Contents: Gene symbol. Order matters - add new genes to bottom")
	Temp=vcat(Temp,"// Records: kegg_organism_code:gene_symbol.splice")
	Temp=vcat(Temp,"// Records: hsa:10195.1")
	Temp=vcat(Temp,"#GENE-SYMBOL-ORDER::START ------------------------------------ //")
	GSO=Gene_Symbol_Order_function(RGMVFF)	
	Temp=vcat(Temp,GSO)
	Temp=vcat(Temp,"#GENE-SYMBOL-ORDER::END -------------------------------------- //")

	#writing .VFF to .txt file
	writedlm(Save_VFF_name,Temp)
	end


  @doc """	
	add_reaction(VFF_name::String,reaction_name::String,reaction_descriptor::String,reaction_representation::String,lowerbound:Float64,upperbound::Float64,decoding_list::String,VFF_new_name::String,KOC::String,rule::String) \n
	VFF_name-> Name of VFF_name file which one wants to add reaction to \n
	reaction_name-> Name of reaction \n
	reaction_descriptor-> Description of reaction \n
	reaction_representation-> The symbolic representation of the reaction \n
	lowerbound-> The lowerbound for the flux for the added reaction \n
	upperbound-> The upperbound for the flux for the added reaction \n
	decoding_list-> The symbolic molecules in the reaction_representation and there corresponding formal name, in the following form 'symbolic name=formal name,symbolic name=formal name, ...' \n
	VFF_new_name-> The name of the .txt file one wants to create and save to, should have .txt in its name \n
	KOC-> What organism are you dealing with and inputed as a Kegg organism code \n
	rule-> The new rule / the rule associated with the reaction one is addeding to the reaction network \n

	This function adds on a reaction, the corresponding metabolites, and gene sections to the correct sections by taking in a previously constructed .VFF.

       """ ->
	function add_reaction(VFF_name::String,reaction_name::String,reaction_descriptor::String,reaction_representation::String,lowerbound::Float64,upperbound::Float64,decoding_list::String,VFF_new_name::String,KOC::String,rule::String)
	    #Opening the file to be edited
            f=open(VFF_name)
	    lines=readlines(f)

	    #seperation of data
	    rxn_names=reaction_name
	    rxns=reaction_representation

	    #building reactions list
	    reactants_list=Array{Union{Nothing, String}}(nothing, 1, 1)
	    products_list=Array{Union{Nothing, String}}(nothing, 1, 1)
	    for i in 1:1#length(rxns) #column i of matrix
		bigsplit=split(rxns,r"  -> |  <=> ")#splits rxn into left and right half
		left_rxn=bigsplit[1]
		left_chem_raw=split(left_rxn," + ")#subdivides the left half into chemicals and their coeffiecent
		left_chem_and_co=split.(left_chem_raw," ")#splits up chemicals and their coeffiecent
		reactants=""
		products=""
		for j in 1:length(left_chem_and_co)
			if j!=length(left_chem_and_co)       
			    if length(left_chem_and_co[j])==2 && left_chem_and_co[j][2]!="" #if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
				coeff=left_chem_and_co[j][1]
				chemical=left_chem_and_co[j][2]			
				reactants=reactants*"$coeff"*"*"*"$chemical"*"+"
			    elseif length(left_chem_and_co[j])==1 && left_chem_and_co[j][1]!=""#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
				chemical=left_chem_and_co[j][1]			
				reactants=reactants*"$chemical"*"+"
			    end
			else
			    if length(left_chem_and_co[j])==2 && left_chem_and_co[j][2]!="" #if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
				coeff=left_chem_and_co[j][1]
				chemical=left_chem_and_co[j][2]			
				reactants=reactants*"$coeff"*"*"*"$chemical"
			    elseif length(left_chem_and_co[j])==1 && left_chem_and_co[j][1]!=""#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
				chemical=left_chem_and_co[j][1]			
				reactants=reactants*"$chemical"
			    end
			end
		end
		reactants_list[i]=reactants
		right_rxn=bigsplit[2]
		right_chem_raw=split(right_rxn," + ")#subdivides the right half into chemicals and their coeffiecent
		right_chem_and_co=split.(right_chem_raw)#splits up chemicals and their coeffiecent
		for j in 1:length(right_chem_and_co)
			if j!=length(right_chem_and_co)          
			   if length(right_chem_and_co[j])==2 && right_chem_and_co[j][2]!=""#if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
				coeff=right_chem_and_co[j][1]
				chemical=right_chem_and_co[j][2]			
				products=products*"$coeff"*"*"*"$chemical+"
			    elseif length(right_chem_and_co[j])==1 && right_chem_and_co[j][1]!=""#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
				chemical=right_chem_and_co[j][1]			
				products=products*"$chemical+"
			    end
			else
			    if length(right_chem_and_co[j])==2 && right_chem_and_co[j][2]!=""#if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
				coeff=right_chem_and_co[j][1]
				chemical=right_chem_and_co[j][2]			
				products=products*"$coeff"*"*"*"$chemical"
			    elseif length(right_chem_and_co[j])==1 && right_chem_and_co[j][1]!=""#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
				chemical=right_chem_and_co[j][1]			
				products=products*"$chemical"
			    end
			end	
		end
		products_list[i]=products
	    end

	    #creation of V_boundaries_vector
	    V_boundaries_vector=[lowerbound upperbound]

	    #Seperating up data
	    Descriptor_list=reaction_descriptor

	    #Initializing E.C. Storage List
	    ec_list=Array{Union{Nothing, String}}(nothing, 1, 1)

	    #obtaining E.C list
	    for v in 1:1
		ec_temp=[]
		Descriptor=Descriptor_list
	    	pos_ecs=read(`curl -g http://rest.kegg.jp/find/enzyme/$Descriptor/`, String)
		sep_pos_ecs=split(pos_ecs,"\n")
			for u in 1:length(sep_pos_ecs)
				s_ec=split(sep_pos_ecs[u],"\t")
			   if length(s_ec)!=1
				r_descriptors=split(s_ec[2],"; ")#retrived discriptors of enzymes
				for w in 1:1
					if isequal(Descriptor,r_descriptors)
	 					ec_temp=push!(ec_temp,s_ec[1])
					end
				end
			   end
			end
		ec_for_rxn=""	
		for s in 1:length(ec_temp)
			temp=ec_temp[s]		
			ec_for_rxn="$ec_for_rxn"*"$temp"*"::"
		end
		ec_for_rxn=ec_for_rxn*"[]"
		ec_list=ec_for_rxn	
	    end
	    
	    #Initializing reaction .vff Storage List
	    rxn_vff_list=Array{Union{Nothing, String}}(nothing, 1, 1)

		
	    #constructing reaction .vff string
	    for t in 1:1
		rxn_temp=rxn_names
		ec_temp=ec_list
		reactant_temp=reactants_list[t]
		product_temp=products_list[t]
		lb_temp=V_boundaries_vector[t,1]
		ub_temp=V_boundaries_vector[t,2]

		#constructing full string
		rxn_vff_temp="$rxn_temp,"*"$ec_temp,"*"$reactant_temp,"*"$product_temp,"*"$lb_temp,"*"$ub_temp"
		
		#saving constructed string
		rxn_vff_list[t]=rxn_vff_temp
	    end

		tt=1
		while tt<length(lines)
			if lines[tt]=="#REACTION::END ----------------------------------------------- //"
				break
			end
		tt=tt+1
		end

		#Updating reaction section
		new_lines=vcat(vcat(lines[1:(tt-1)],rxn_vff_list),lines[tt:length(lines)])

	    #constructing the new metabolite VFF lines
	    meta=split(decoding_list,",")

	    #Prebuilding (symbolic name, real name, kegg ID)
	    new_VFF_meta=Array{Union{Nothing, String}}(nothing, length(meta), 3)
	    
	    for g in 1:length(meta)
		#sorting of chemicals by their compartmental location
		#decoding the symboles and then quering kegg
		new_VFF_meta[g,1]=split(meta[g],"=")[1]
		new_VFF_meta[g,2]="nothing"
		new_VFF_meta[g,3]="nothing"
		        if occursin("[c]",split(meta[g],"=")[1])
		            real_name=split(split(meta[g],"=")[2]," cytosol")[1]
		            new_VFF_meta[g,2]=real_name
			    real_name=replace(real_name, " " => "%20") ###correction from previous version
		            h=read(`curl -g http://rest.kegg.jp/find/compound/$real_name/`, String)
		            name=split(h,"\t")
				if length(name)>=2
				    Cid=name[1]
				    new_VFF_meta[g,3]=Cid
				end
		        elseif occursin("[m]",split(meta[g],"=")[1])
		            real_name=split(split(meta[g],"=")[2]," mitochondria")[1]
		            new_VFF_meta[g,2]=real_name
			    real_name=replace(real_name, " " => "%20") ###correction from previous version
		            h=read(`curl -g http://rest.kegg.jp/find/compound/$real_name/`, String)
		            name=split(h,"\t")
				if length(name)>=2
				    Cid=name[1]
				    new_VFF_meta[g,3]=Cid
				end
		        elseif occursin("[r]",split(meta[g],"=")[1])
		            real_name=split(split(meta[g],"=")[2]," ER")[1]
		            new_VFF_meta[g,2]=real_name
			    real_name=replace(real_name, " " => "%20") ###correction from previous version
		            h=read(`curl -g http://rest.kegg.jp/find/compound/$real_name/`, String)
		            name=split(h,"\t")
				if length(name)>=2
				    Cid=name[1]
				    new_VFF_meta[g,3]=Cid
				end
		        elseif occursin("[x]",split(meta[g],"=")[1])
		            real_name=split(split(meta[g],"=")[2]," peroxisome")[1]
		            new_VFF_meta[g,2]=real_name
			    real_name=replace(real_name, " " => "%20") ###correction from previous version
		            h=read(`curl -g http://rest.kegg.jp/find/compound/$real_name/`, String)
		            name=split(h,"\t")
				if length(name)>=2
				    Cid=name[1]
				    new_VFF_meta[g,3]=Cid
				end
		        elseif occursin("[l]",split(meta[g],"=")[1])
		            real_name=split(split(meta[g],"=")[2]," lysome")[1]
		            new_VFF_meta[g,2]=real_name
			    real_name=replace(real_name, " " => "%20") ###correction from previous version
		            h=read(`curl -g http://rest.kegg.jp/find/compound/$real_name/`, String)
		            name=split(h,"\t")
				if length(name)>=2
				    Cid=name[1]
				    new_VFF_meta[g,3]=Cid
				end
		        elseif occursin("[e]",split(meta[g],"=")[1])
		            real_name=split(split(meta[g],"=")[2]," extracellular")[1]
		            new_VFF_meta[g,2]=real_name
			    real_name=replace(real_name, " " => "%20") ###correction from previous version
		            h=read(`curl -g http://rest.kegg.jp/find/compound/$real_name/`, String)
		            name=split(h,"\t")
				if length(name)>=2
				    Cid=name[1]
				    new_VFF_meta[g,3]=Cid
				end
		        end
	    end

		#Pre-construting new VFF array to be saved to external file which can be later read and then constructed into the complete .VFF for the corresponding model
		VFF_meta_array=Array{Union{Nothing, String}}(nothing, length(new_VFF_meta[:,1]), 1)

		for b in 1:length(VFF_meta_array[:,1])

			#assigning symbolic name
			sym_name=new_VFF_meta[b,1]
			
			#assigning real name
			real_name=new_VFF_meta[b,2]	

			#assigning Kegg ID
			Kegg_ID=new_VFF_meta[b,3]

			#Concatenating into correct .VFF format
			VFF_metabolite="$sym_name"*"="*"$real_name"*"::"*"$Kegg_ID"

			#Adding it to the VFF_metabolite_array
			VFF_meta_array[b]=VFF_metabolite
		end

	    #Finding metabolites section
	    tt=1
	    while tt<length(new_lines)
			if new_lines[tt]=="#METABOLITES::START ------------------------------------------ //"
				break
			end
		tt=tt+1
	    end
		
	    jj=1
	    while jj<length(new_lines)
			if new_lines[jj]=="#METABOLITES::END -------------------------------------------- //"
				break
			end
		jj=jj+1
	    end

	    #Grabbing the metabolite section
	    v=new_lines[tt+1:jj-1]
		
	    #Updating the metabolite section
	    total_meta=union(v,VFF_meta_array)

	    #prebuilding of chemical's compartment organization vectors used to organize the total_meta section
	    c_unsort=[]
	    m_unsort=[]
	    r_unsort=[]
	    x_unsort=[]
	    l_unsort=[]
	    e_unsort=[]

	    for i in 1:length(total_meta)#length(rxns) #column i of matrix
		#Bin sorting of chemicals by their compartmental location
		if occursin("[c]",total_meta[i])
		    c_unsort=push!(c_unsort,total_meta[i])
		elseif occursin("[m]",total_meta[i])
		    m_unsort=push!(m_unsort,total_meta[i])
		elseif occursin("[r]",total_meta[i])
		    r_unsort=push!(r_unsort,total_meta[i])
		elseif occursin("[x]",total_meta[i])
		    x_unsort=push!(x_unsort,total_meta[i])
		elseif occursin("[l]",total_meta[i])
		    l_unsort=push!(l_unsort,total_meta[i])
		elseif occursin("[e]",total_meta[i])
		    e_unsort=push!(e_unsort,total_meta[i])
		end
	    end

	    #Sorting chemicals via ASCII assignments
	    c_sorted=sort(c_unsort)
	    m_sorted=sort(m_unsort)
	    r_sorted=sort(r_unsort)
	    x_sorted=sort(x_unsort)
	    l_sorted=sort(l_unsort)
	    e_sorted=sort(e_unsort)

	    #Prebuilding sorted chemical vector
	    chemicals_sorted=[]

	    #Building sorted chemical vector by adding blocks of sorted chemicals
	    chemicals_sorted=c_sorted
	    chemicals_sorted=append!(chemicals_sorted,m_sorted)
	    chemicals_sorted=append!(chemicals_sorted,r_sorted)
	    chemicals_sorted=append!(chemicals_sorted,x_sorted)
	    chemicals_sorted=append!(chemicals_sorted,l_sorted)
	    chemicals_sorted=append!(chemicals_sorted,e_sorted)

	    #Updating VFF
	    VFF_updated=vcat(vcat(new_lines[1:tt],chemicals_sorted),new_lines[jj:length(new_lines)])

	    #Updating Gene sections
		#Constructing the information and formating it for the the array then storing it
		RGMVFF=""
		for i in 1:1
		name=reaction_name
		raw_data=read(`curl -g http://rest.kegg.jp/find/genes/$name/`, String)
		s_data=split(raw_data,"\n")
		gene_storage=""
			for j in 1:length(s_data)
			ss_data=split(s_data[j],"\t")
				if length(ss_data)>1
					if occursin(name,ss_data[2])
						if occursin(KOC,ss_data[1])
							gene_storage=gene_storage*ss_data[1]*".1"*"::" #Assuming normal unspecial spliced gene if one wants special splices one has to hand edit the gene section
						end
					end
				end 
			end
		RGMVFF=name*"="*gene_storage*"[]"
		end

	    #Finding information needed to determine lines to which append
	    gg=1
	    while gg<length(VFF_updated)
			if VFF_updated[gg]=="#REACTION-GENE-MAP::START ------------------------------------ //"
				break
			end
		gg=gg+1
	    end
		
	    hh=1
	    while hh<length(VFF_updated)
			if VFF_updated[hh]=="#REACTION-GENE-MAP::END -------------------------------------- //"
				break
			end
		hh=hh+1
	    end

	    #Grabbing and updating reaction gene map section
	    v=vcat(VFF_updated[gg+1:hh-1],RGMVFF)

	    #Updating VFF
	    New_VFF_updated=vcat(vcat(VFF_updated[1:gg],v),VFF_updated[hh:length(VFF_updated)])

	    #Updating Gene-Symbol-section
		#Preconstructing Gene_symbol_order .VFF array
		GSO=[]
		
		#Parsing and then building GSO
		for i in 1:1
			Genes=split(RGMVFF,"=")
			I_Genes=split(Genes[2],"::")
			for j in 1:length(I_Genes)
				if I_Genes[j]!="[]"
					GSO=push!(GSO,I_Genes[j])
				end
			end
		end

	    #Finding information needed to determine line which to append
	    xx=1
	    while xx<length(VFF_updated)
			if VFF_updated[xx]=="#GENE-SYMBOL-ORDER::START ------------------------------------ //"
				break
			end
		xx=xx+1
	    end
		
	    zz=1
	    while zz<length(VFF_updated)
			if VFF_updated[zz]=="#GENE-SYMBOL-ORDER::END -------------------------------------- //"
				break
			end
		zz=zz+1
	    end

	    #Updating Gene symbol order section
	    v=vcat(New_VFF_updated[xx+1:zz-1],GSO)

	    #Updating VFF
	    New_VFF_updated=vcat(vcat(New_VFF_updated[1:xx],v),New_VFF_updated[zz:length(New_VFF_updated)])

		#Finding rules section
		a=1
		while a<length(New_VFF_updated)
			if New_VFF_updated[a]=="#RULES::START ------------------------------------------------ //"
				break
			end
		a=a+1
		end

		b=1
		while b<length(New_VFF_updated)
			if New_VFF_updated[b]=="#RULES::END -------------------------------------------------- //"
				break
			end
		b=b+1
		end

		#Grabbing rules section
		rule_data=New_VFF_updated[a+1:b-1]

		#Creating the new rule line
		new_rule_line=reaction_name*"="*rule

		#Updating the rule data
		v=push!(rule_data,new_rule_line)

		#Updating .VFF
		New_VFF_updated=vcat(vcat(New_VFF_updated[1:a],v),New_VFF_updated[b:length(New_VFF_updated)])

	    #writing .VFF to .txt file
	    writedlm(VFF_new_name,New_VFF_updated)#Save new file
	end

  @doc """	
	delete_reaction(VFF_name::String,reaction_name::String,VFF_new_name::String) \n
	VFF_name-> Name of VFF_name file which one wants to remove a reaction from \n
	reaction_name-> Name of reaction that one wants to remove \n
	VFF_new_name-> The name of the .txt file one wants to create and save to, should have .txt in its name \n

	This function removes a reaction and corrects the corresponding metabolites sections by taking in a previously constructed .VFF and the reaction one wants to remove.

       """ ->
	function delete_reaction(VFF_name::String,reaction_name::String,VFF_new_name::String)
		#Opening the file to be edited	
		f=open(VFF_name)
		lines=readlines(f)

		#Finding the location of the deleted reaction in the .VFF
		t=0
		for i in 1:length(lines)
			if isequal(reaction_name,split(lines[i],",")[1])
				t=i
				break
			end
		end

		#Saving the deleted reactions information for use in eliminating now obsolet metabolites and such
		del_rxn=lines[t]

		#Rebuilding the .VFF without the old reaction
		lines=vcat(lines[1:t-1],lines[t+1:length(lines)])

		#Eliminating the old metabolites that are no longer present in the reaction network
		pos_met=split(del_rxn,",")[3:4]

		#Intializing an array to store the metabolites that are possible on the chopping block
		store=[]
		
		#Creating list of the metabolites on the chopping block
		for j in 1:length(pos_met)
			met=split(pos_met[j],"+")
				for k in 1:length(met)
					meta_temp=split(met[k],"*")
						if length(meta_temp)==2
							store=push!(store,meta_temp[2])
						else
							store=push!(store,meta_temp[1])
						end
				end
		end

		#Finding the reaction section
		tt=1
		while tt<length(lines)
			if lines[tt]=="#REACTION::START --------------------------------------------- //"
				break
			end
		tt=tt+1
		end
		
		jj=1
		while jj<length(lines)
			if lines[jj]=="#REACTION::END ----------------------------------------------- //"
				break
			end
		jj=jj+1
		end

		#Finding the metabolite section
		ss=1
		while ss<length(lines)
			if lines[ss]=="#METABOLITES::START ------------------------------------------ //"
				break
			end
		ss=ss+1
		end

		ii=1
		while ii<length(lines)
			if lines[ii]=="#METABOLITES::END -------------------------------------------- //"
				break
			end
		ii=ii+1
		end

		#Building result's storage vector
		results=zeros(length(store),1)

		#Determining which metabolites get removed	
		zz=0
		for l in tt+1:jj-1
			rem_met=split(lines[l],",")[3:4]
			for m in 1:length(rem_met)
				check_met=split(rem_met[m],"+")
				for n in 1:length(check_met)
					metabolite=split(check_met[n],"*")
					if length(metabolite)==2
						for o in 1:length(store)
							if isequal(store[o],metabolite[2])
								results[o]=results[o]+1
							end
						end
					else
						for o in 1:length(store)
							if isequal(store[o],metabolite[1])
								results[o]=results[o]+1
							end
						end
					end

				end
			end
		end

		#Setting up removal vector
		removal=[]

		#Determining the lines to be removed
		for p in 1:length(results)
			if results[p]==0
				for q in ss+1:ii-1
					tested_met=split(lines[q],"=")[1]
					if isequal(store[p],tested_met)
						removal=append!(removal,q)
						break
					end	
				end
			end
		end

		#Sorting removal so that the metabolites can be removed in an easy accurate manner	
		org_removal=sort(removal)
		
		#Removing the outdated metabolites
		for r in 1:length(org_removal)
			s=length(org_removal)-r+1
			lines=vcat(lines[1:org_removal[s]-1],lines[org_removal[s]+1:length(lines)])
		end

		#Updating gene section

		#Finding information needed to determine lines to which append
		gg=1
		while gg<length(lines)
			if lines[gg]=="#REACTION-GENE-MAP::START ------------------------------------ //"
				break
			end
		gg=gg+1
		end

		hh=1
		while hh<length(lines)
			if lines[hh]=="#REACTION-GENE-MAP::END -------------------------------------- //"
				break
			end
		hh=hh+1
		end

		#Building the curation vector
		vold=lines[gg+1:hh-1]
		
		#Finding the location of the deleted reaction's gene mapping in the .VFF
		t=0
		for i in gg:length(lines)
			if isequal(reaction_name,split(lines[i],"=")[1])
				t=i
				break
			end
		end
		
		#Removing the bad genes
		vnew=vcat(lines[gg+1:t-1],lines[t+1:hh-1])
		
		#Updating reaction-gene-section
		lines=vcat(vcat(lines[1:gg+1],vnew),lines[hh:length(lines)])

		#Preconstructing Gene_symbol_order .VFF array
		GSO_new=[]
		
		#Parsing and then building GSO
		for i in 1:length(vnew)
			Genes=split(vnew[i],"=")
			I_Genes=split(Genes[2],"::")
			for j in 1:length(I_Genes)
				if I_Genes[j]!="[]"
					GSO_new=push!(GSO_new,I_Genes[j])
				end
			end
		end

		#Finding information needed to determine line with which one can use to correct the .VFF and specifically the Gene Symbolic Order
		xx=1
		while xx<length(lines)
			if lines[xx]=="#GENE-SYMBOL-ORDER::START ------------------------------------ //"
				break
			end
		xx=xx+1
		end

		zz=1
		while zz<length(lines)
			if lines[zz]=="#GENE-SYMBOL-ORDER::END -------------------------------------- //"
				break
			end
		zz=zz+1
		end

		#Updating .VFF
		lines=vcat(vcat(lines[1:xx],GSO_new),lines[zz:length(lines)])

		#Finding rules section
		a=1
		while a<length(lines)
			if lines[a]=="#RULES::START ------------------------------------------------ //"
				break
			end
		a=a+1
		end

		b=1
		while b<length(lines)
			if lines[b]=="#RULES::END -------------------------------------------------- //"
				break
			end
		b=b+1
		end

		#Grabbing rules section
		rule_data=lines[a+1:b-1]

		#Finding the location of the deleted reaction in the Rule section in the .VFF
		t=0
		for i in 1:length(rule_data)
			if isequal(reaction_name,split(rule_data[i],"=")[1])
				t=i
				break
			end
		end

		#Removing rule associated with the removed reaction
		v=vcat(rule_data[1:t-1],rule_data[t+1:length(rule_data)])

		#Initializing deletion reaction array
		del_array=[]

		#As the complexity of the rules section can very significantly it is left to the user to edit the rules section but this function informs the user of the lines in which must be
	 	#removed by hand, the reaction number is replaced by and *
		#The following block corrects the change in reaction number		
		for i in 1:length(v)
			for j in 1:length(v)+1
				temp=split(v[i],"($j)")	
					if j==t
						if length(temp)==2
							v[i]=temp[1]*"("*"*"*")"*temp[2]
							del_index=i+a							
							del_array=push!(del_array,del_index)
						end
					elseif j>t
						if length(temp)==2
							jj=j-1
							v[i]=temp[1]*"("*"$jj"*")"*temp[2]
						end
					end
			end
		end

		#Telling the user which lines must have their logic updated
		println("The following lines are in the rules section and the outputed line numbers are refering to the total line number, and these lines need to have their logic updated and are printed in the following array:")

		#Telling the user the lines
		println(del_array)

		#Updating the .VFF lines
		lines=vcat(vcat(lines[1:a],v),lines[b:length(lines)])

		#Saving the updated .VFF	
		writedlm(VFF_new_name,lines)#Save new file
	end


  @doc """	
	adjust_reaction(VFF_name::String,reaction_name::String,New_lowerbound::Float64,New_upperbound::Float64,VFF_new_name::String) \n
	VFF_name-> Name of VFF_name file which one wants to adjust the flux of a reaction from \n
	reaction_name-> Name of reaction that one wants to adjust \n
	New_lowerbound-> The new lowerbound for the flux of the reaction of interest \n
	New_upperbound-> The new upperbound for the flux of the reaction of interest \n
	VFF_new_name-> The name of the .txt file one wants to create and save to, should have .txt in its name \n

	This function adjusts the flux's of an reaction of interest in the .VFF, then saves the adjusted .VFF.

       """ ->
	function adjust_reaction(VFF_name::String,reaction_name::String,New_lowerbound::Float64,New_upperbound::Float64,VFF_new_name::String)
		#Opening the file to be edited	
		f=open(VFF_name)
		lines=readlines(f)

		#Finding the location of the reaction to be adjust in the .VFF
		t=0
		for i in 1:length(lines)
			if isequal(reaction_name,split(lines[i],",")[1])
				t=i
				break
			end
		end
		
		#Deliminating reaction to allow for easy editing	
		broken_up_line=split(lines[t],",")
		
		#Adjusting the reaction
		broken_up_line[5]="$New_lowerbound"
		broken_up_line[6]="$New_upperbound"

		#Reconstructing full string and replacing old reaction string with adjusted one
		lines[t]="$(broken_up_line[1]),"*"$(broken_up_line[2]),"*"$(broken_up_line[3]),"*"$(broken_up_line[4]),"*"$(broken_up_line[5]),"*"$(broken_up_line[6])"

		#Saving the newly adjust .VFF
		writedlm(VFF_new_name,lines)#Save new file
	end


  @doc """	
	Reaction_Gene_Map_function(File_Name_txt::String,KOC::String,first_rxn_row::Int,last_rxn_row::Int,rxn_ref_name_column::Int) \n
	File_name_txt-> Name of the file containing the reactions and the corresponding information \n
	KOC-> What organism are you dealing with and inputed as a Kegg organism code \n
	first_rxn_row-> The first row containing a reaction in the .txt file \n 
	last_rxn_row-> The last row containing a reaction in the .txt file \n
	rxn_ref_name_column-> The column number which contains the reaction name column \n

	This function adjusts the flux's of an reaction of interest in the .VFF, then saves the adjusted .VFF.

       """ ->
	function Reaction_Gene_Map_function(File_Name_txt::String,KOC::String,first_rxn_row::Int,last_rxn_row::Int,rxn_ref_name_column::Int)

		#Collection of data
		rxn_path_data=readdlm(File_Name_txt,'\t',String,'\n')

		#Seperation of data
		rxn_names=rxn_path_data[first_rxn_row:last_rxn_row,rxn_ref_name_column]#finds how rxn are named used for reference

		#Preconstructing storage vector to store the constructed .VFF lines for the Reaction-Gene-Mapping section
		RGMVFF=Array{Union{Nothing, String}}(nothing, length(rxn_names), 1)

		#Constructing the information and formating it for the the array then storing it
		for i in 1:length(rxn_names)
		name=rxn_names[i]
		raw_data=read(`curl -g http://rest.kegg.jp/find/genes/$name/`, String)
		s_data=split(raw_data,"\n")
		gene_storage=""
			for j in 1:length(s_data)
			ss_data=split(s_data[j],"\t")
				if length(ss_data)>1
					if occursin(name,ss_data[2])
						if occursin(KOC,ss_data[1])
							gene_storage=gene_storage*ss_data[1]*".1"*"::" #Assuming normal unspliced gene if spliced gene is added one must hand adjust it
						end
					end
				end 
			end
		RGMVFF[i]=name*"="*gene_storage*"[]"
		end

		#Returning a reaction gene mapping vector in .VFF format		
		return RGMVFF
	end


  @doc """	
	Gene_Symbol_Order_function(RGMVFF::Array{Union{Nothing, String},2})\n
	RGMVFF-> The string array containing the reaction gene mapping information\n

	This function creates the Gene_Symbol_Order section string array for the .VFF construction.

       """ ->
	function Gene_Symbol_Order_function(RGMVFF::Array{Union{Nothing, String},2})
		#Preconstructing Gene_symbol_order .VFF array
		GSO=[]
		
		#Parsing and then building GSO
		for i in 1:length(RGMVFF)
			Genes=split(RGMVFF[i],"=")
			I_Genes=split(Genes[2],"::")
			for j in 1:length(I_Genes)
				if I_Genes[j]!="[]"
					GSO=push!(GSO,I_Genes[j])
				end
			end
		end
		
		#Returning a gene symbol order vector in .VFF format
		return GSO
	end

  @doc """	
	Rule_building(File_Name_txt::String,rxns_column_number::Int,first_rxn_row::Int,last_rxn_row::Int,rxn_ref_name_column::Int,Rule_Document::String,Rule_String_Array::Array) \n
	File_Name_txt-> Where the reactions are saved in a tab deliminated .txt file \n
	rxns_column_number-> The column in which the symbolic representation of the reaction is located \n
	first_rxn_row-> The first row containing a reaction in the .txt file \n
	last_rxn_row-> The last row containing a reaction in the .txt file \n
	rxn_ref_name_column-> The column number which contains the reaction name column \n
	Rule_Document-> The document location / name which contians the information of the rules for the reactions \n
	Rule_String_Array-> A String array that contains the rule information in .VFF format \n

	This function creates / builds a rule string array that is in .VFF format.

       """ ->
	function Rule_building(File_Name_txt::String,rxns_column_number::Int,first_rxn_row::Int,last_rxn_row::Int,rxn_ref_name_column::Int,Rule_Document::String,Rule_String_Array::Array)
		
		#Checking if Rule constraint are given in either of the two acceptable forms i.e a string array that just needs to be append to the .VFF 
		if Rule_Document!=""
			#Collection of data
			Rule_Vector=readdlm(Rule_Document, '\t', String, '\n')[:]
			return Rule_Vector
		elseif Rule_String_Array!=[]
			return Rule_String_Array
		else
		#Now having to construct and empty rule section for later
			#Pre-intializing Rule Section
			Rule_Section=Array{String,1}()

			#Collection of data
			rxn_path_data=readdlm(File_Name_txt, '\t', String, '\n')

			#Seperation of data
			rxn_names=rxn_path_data[first_rxn_row:last_rxn_row,rxn_ref_name_column]#finds how rxn are named used for reference	
			
			#Building empty Rule Section
			for i in 1:length(rxn_names)
				temp=rxn_names[i]*"="*"[]"
				Rule_Section=push!(Rule_Section,temp)
			end
			
			#Returning the newly created Rule Section
			return Rule_Section
			
		end

	end

  @doc """	
	Rule_edit(VFF_name::String,reaction_name::String,New_rule::String,New_VFF_name::String) \n
	VFF_name-> This is the .VFF which one wants the edit a rule in \n
	reaction_name-> This is the name of the reaction which one wants to edit the rule of the reaction \n
	New_rule-> This is the new rule for the requested reaction \n
	New_VFF_name-> This is the name which one wants the new file saved under \n

	This function edits the rules associated with a desired reaction.

       """ ->
	function Rule_edit(VFF_name::String,reaction_name::String,New_rule::String,New_VFF_name::String)
		#Opening the file to be edited
		f=open(VFF_name)
		lines=readlines(f)

		#Finding information needed to determine line with which one can use to correct the .VFF and specifically the Rules
		a=1
		while a<length(lines)
			if lines[a]=="#RULES::START ------------------------------------------------ //"
				break
			end
		a=a+1
		end

		b=1
		while b<length(lines)
			if lines[b]=="#RULES::END -------------------------------------------------- //"
				break
			end
		b=b+1
		end

		#Creating editing array
		v=lines[a+1:b-1]

		#Searching for the correct reaction		
		t=1
		while t<=length(v)
			if split(v[t],"=")[1]==reaction_name
				new_line=reaction_name*"="*New_rule
				v[t]=new_line
				break			
			end
			t=t+1
		end
			
		#Correcting the document
		new_lines=vcat(vcat(lines[1:a],v),lines[b:length(lines)])		

		#Saving the newly adjust .VFF
		writedlm(New_VFF_name,new_lines)#Save new file	
	end

  @doc """	
	VFF_to_mat(VFF_file_name::String,New_file_name::String,model_name::String,constraint::Array,File_Name_txt::String,rxns_column_number::Int,first_rxn_row::Int,last_rxn_row::Int,rxn_ref_name_column::Int,rxn_des_column::Int,subsystems_column::Int,b::Array) \n
	VFF_file_name-> The file name of the .VFF that one wants to convert into the .mat file \n
	New_file_name-> The name that one wants to assigned to the new .mat file \n
	model_name-> The name that one wants the containing structure in MATLAB \n
	constraint-> The constraint c that one is optimizing the network with respect to \n
	File_Name_txt-> The file that contians the reaction network from which one is going to use to build the stoichiometric matrix \n
	rxns_column_number-> The column in which the symbolic representation of the reaction is located \n
	first_rxn_row-> The first row containing a reaction in the .txt file \n
	last_rxn_row-> The last row containing a reaction in the .txt file \n
	rxn_ref_name_column-> The column number which contains the reaction name column \n
	rxn_des_column-> The column number which contains the description information for each reaction \n 
	subsystems_column-> The column number which contains the subsystem information \n
	b-> The constraint vector for which Sv must be equal to \n
	
	This function converts a .VFF for an arbitary cellular system model into the .mat version which is compatable with the cobra toolbox

       """ ->
	function VFF_to_mat(VFF_file_name::String,New_file_name::String,model_name::String,constraint::Array,File_Name_txt::String,rxns_column_number::Int,first_rxn_row::Int,last_rxn_row::Int,rxn_ref_name_column::Int,rxn_des_column::Int,subsystems_column::Int,b::Array)

		lines=readlines(open(VFF_file_name))

		#Finding the reaction section
		tt=1
		while tt<length(lines)
			if lines[tt]=="#REACTION::START --------------------------------------------- //"
				break
			end
		tt=tt+1
		end

		jj=1
		while jj<length(lines)
			if lines[jj]=="#REACTION::END ----------------------------------------------- //"
				break
			end
		jj=jj+1
		end

		#Grabbing reaction section
		rxn_sec=lines[tt+1:jj-1]
		
		#Preconstructing vectors
		rxns_m=Array{String,1}()	
		mets=Array{String,1}()
		metNames=Array{String,1}()
		rules=Array{String,1}()
		rev=Array{Float64,1}()
		lb=Array{Float64,1}()
		ub=Array{Float64,1}()
		genes=Array{String,1}()


		#Constructing rxns, rev, lb, ub vectors
		for i in 1:length(rxn_sec)
			rxns_m=push!(rxns_m,split(rxn_sec[i],",")[1])
				if parse(Float64,split(rxn_sec[i],",")[5])<0
					if parse(Float64,split(rxn_sec[i],",")[6])>0
						rev=push!(rev,1)
					else
						rev=push!(rev,0)
					end
				else
					rev=push!(rev,0)
				end
			lb=push!(lb,parse(Float64,split(rxn_sec[i],",")[5]))
			ub=push!(ub,parse(Float64,split(rxn_sec[i],",")[6]))	
		end

		#Building the arbitary bias optimization constraint if it is not specified
		if length(constraint)==0
			c=zeros(length(rxn_sec),1)
		else
			c=constraint	
		end

		#Finding the gene section
		xx=1
		while xx<length(lines)
			if lines[xx]=="#GENE-SYMBOL-ORDER::START ------------------------------------ //"
				break
			end
		xx=xx+1
		end

		zz=1
		while zz<length(lines)
			if lines[zz]=="#GENE-SYMBOL-ORDER::END -------------------------------------- //"
				break
			end
		zz=zz+1
		end

		#Grabbing gene section
		gene_sec=lines[xx+1:zz-1]

		#Building gene section
		genes=append!(genes,gene_sec)

		#Finding the metabolite section
		ss=1
		while ss<length(lines)
			if lines[ss]=="#METABOLITES::START ------------------------------------------ //"
				break
			end
		ss=ss+1
		end

		ii=1
		while ii<length(lines)
			if lines[ii]=="#METABOLITES::END -------------------------------------------- //"
				break
			end
		ii=ii+1
		end

		#Grabbing the metabolite data
		meta_sec=lines[ss+1:ii-1]

		#Building the arbitary constraint if it is not specified
		if length(b)==0
			b=zeros(length(meta_sec),1)
		else
			b=b
		end

		#Building metabolite section
		for i in 1:length(meta_sec)
			mets=push!(mets,split(meta_sec[i],"=")[1])
			metNames=push!(metNames,split(split(meta_sec[i],"=")[2],"::")[1])
		end

		#Pre-building reaction-gene-mapping matrix
		rxnGeneMat=zeros(length(rxn_sec),length(gene_sec))

		#Finding RGM section
		gg=1
		while gg<length(lines)
			if lines[gg]=="#REACTION-GENE-MAP::START ------------------------------------ //"
				break
			end
		gg=gg+1
		end

		hh=1
		while hh<length(lines)
			if lines[hh]=="#REACTION-GENE-MAP::END -------------------------------------- //"
				break
			end
		hh=hh+1
		end

		#Obtaining rgm section
		rgm_sec=lines[gg+1:hh-1]

		#Building rxnGeneMat
		for i in 1:length(rgm_sec)
			for j in 1:length(gene_sec)
				if occursin(gene_sec[j]*"::",rgm_sec[i])
					rxnGeneMat[i,j]=1
				end	
			end
		end

		#Converting rxnGeneMat from Float64 to Int64
		rxnGeneMat=convert(Array{Int64},rxnGeneMat)

		#Building S_matrix

		    #Collection of data
		    rxn_path_data=readdlm(File_Name_txt, '\t', String, '\n')

		    #seperation of data
		    rxn_names=rxn_path_data[first_rxn_row:last_rxn_row,rxn_ref_name_column]#finds how rxn are named used for reference
		    rxns=rxn_path_data[first_rxn_row:last_rxn_row,rxns_column_number]#the actual stoichometric reactions
		    subsystems=rxn_path_data[first_rxn_row:last_rxn_row,subsystems_column]
		    rxn_des=rxn_path_data[first_rxn_row:last_rxn_row,rxn_des_column]

		    #building unique metabolic list
		    chemicals_raw=[]
		    for i in 1:length(rxns)#length(rxns) #column i of matrix
			bigsplit=split(rxns[i],r"  -> |  <=> ")#splits rxn into left and right half
			left_rxn=bigsplit[1]
			left_chem_raw=split(left_rxn," + ")#subdivides the left half into chemicals and their coeffiecent
			left_chem_and_co=split.(left_chem_raw," ")#splits up chemicals and their coeffiecent
			for j in 1:length(left_chem_and_co)#collects all unique chemicals on left side
			    if length(left_chem_and_co[j])==2 && (!(in(left_chem_and_co[j][2],chemicals_raw)) && left_chem_and_co[j][2]!="")#if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
				chemicals_raw=push!(chemicals_raw,left_chem_and_co[j][2])
			    elseif length(left_chem_and_co[j])==1 && (!(in(left_chem_and_co[j][1],chemicals_raw)) && left_chem_and_co[j][1]!="")#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
				chemicals_raw=push!(chemicals_raw,left_chem_and_co[j][1])
			    end
			end
			right_rxn=bigsplit[2]
			right_chem_raw=split(right_rxn," + ")#subdivides the right half into chemicals and their coeffiecent
			right_chem_and_co=split.(right_chem_raw)#splits up chemicals and their coeffiecent
			for j in 1:length(right_chem_and_co)#collects all unique chemicals on right side
			    if length(right_chem_and_co[j])==2 && !(in(right_chem_and_co[j][2],chemicals_raw)) && !(right_chem_and_co[j][2]=="")#if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
				chemicals_raw=push!(chemicals_raw,right_chem_and_co[j][2])
			    elseif length(right_chem_and_co[j])==1 && !(in(right_chem_and_co[j][1],chemicals_raw)) && !(right_chem_and_co[j][1]=="")#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
				chemicals_raw=push!(chemicals_raw,right_chem_and_co[j][1])
			    end
			end
		    end

		    #prebuilding of chemical's compartment organization vectors
		    c_unsort=[]
		    m_unsort=[]
		    r_unsort=[]
		    x_unsort=[]
		    l_unsort=[]
		    e_unsort=[]

		    for i in 1:length(chemicals_raw)#length(rxns) #column i of matrix
			#Bin sorting of chemicals by their compartmental location
			if occursin("[c]",chemicals_raw[i])
			    c_unsort=push!(c_unsort,chemicals_raw[i])
			elseif occursin("[m]",chemicals_raw[i])
			    m_unsort=push!(m_unsort,chemicals_raw[i])
			elseif occursin("[r]",chemicals_raw[i])
			    r_unsort=push!(r_unsort,chemicals_raw[i])
			elseif occursin("[x]",chemicals_raw[i])
			    x_unsort=push!(x_unsort,chemicals_raw[i])
			elseif occursin("[l]",chemicals_raw[i])
			    l_unsort=push!(l_unsort,chemicals_raw[i])
			elseif occursin("[e]",chemicals_raw[i])
			    e_unsort=push!(e_unsort,chemicals_raw[i])
			end
		    end

		    #Sorting chemicals via ASCII assignments
		    c_sorted=sort(c_unsort)
		    m_sorted=sort(m_unsort)
		    r_sorted=sort(r_unsort)
		    x_sorted=sort(x_unsort)
		    l_sorted=sort(l_unsort)
		    e_sorted=sort(e_unsort)

		    #Prebuilding sorted chemical vector
		    chemicals_sorted=[]

		    #Building sorted chemical vector by adding blocks of sorted chemicals
		    chemicals_sorted=c_sorted
		    chemicals_sorted=append!(chemicals_sorted,m_sorted)
		    chemicals_sorted=append!(chemicals_sorted,r_sorted)
		    chemicals_sorted=append!(chemicals_sorted,x_sorted)
		    chemicals_sorted=append!(chemicals_sorted,l_sorted)
		    chemicals_sorted=append!(chemicals_sorted,e_sorted)

		    #Creation of base of Stiochiometric matrix
		    S_matrix=zeros(length(chemicals_sorted),length(rxns))

		    #processing of rxns into stoichometric matrix
		    for j in 1:length(rxns)#length(rxns) #column j of matrix
			bigsplit=split(rxns[j],r"  -> |  <=> ")#splits rxn into left and right half
			left_rxn=bigsplit[1]
			left_chem_raw=split(left_rxn," + ")#subdivides the left half into chemicals and their coeffiecent
			left_chem_and_co=split.(left_chem_raw," ")#splits up chemicals and their coeffiecent
			for k in 1:length(left_chem_and_co)#Places coeffiecents in s_matrix and since its the left side the entries recieve an negative sign
			    if length(left_chem_and_co[k])==2 && left_chem_and_co[k][2]!=""#grabs coeffiecent of chemicals if present
				for i in 1:length(chemicals_sorted) #possible row number
				    if occursin(left_chem_and_co[k][2],chemicals_sorted[i])#finds correct row number
				        S_matrix[i,j]=-1*parse(Int,left_chem_and_co[k][1])
				    end
				end
			    elseif length(left_chem_and_co[k])==1 && left_chem_and_co[k][1]!=""#grabs present chemicals and assigns -1 to their S_matrix possition
				for i in 1:length(chemicals_sorted) #possible row number
				    if occursin(left_chem_and_co[k][1],chemicals_sorted[i])#finds correct row number
				        S_matrix[i,j]=-1
				    end
				end
			    end
			end
			right_rxn=bigsplit[2]
			right_chem_raw=split(right_rxn," + ")#subdivides the right half into chemicals and their coeffiecent
			right_chem_and_co=split.(right_chem_raw)#splits up chemicals and their coeffiecent
			for k in 1:length(right_chem_and_co)
			    if length(right_chem_and_co[k])==2 && !(right_chem_and_co[k][2]=="")#grabs coeffiecent of chemicals if present
				for i in 1:length(chemicals_sorted) #possible row number
				    if occursin(right_chem_and_co[k][2],chemicals_sorted[i])#finds correct row number
				        S_matrix[i,j]=parse(Int,right_chem_and_co[k][1])
				    end
				end
			    elseif length(right_chem_and_co[k])==1 && right_chem_and_co[k][1]!=""#grabs present chemicals and assigns 1 to their S_matrix position
				for i in 1:length(chemicals_sorted) #possible row number
				    if occursin(right_chem_and_co[k][1],chemicals_sorted[i])#finds correct row number
				        S_matrix[i,j]=1
				    end
				end
			    end
			end
		    end

		#Renaming S_matrix
		S=S_matrix

		#Finding rules section
		a=1
		while a<length(lines)
			if lines[a]=="#RULES::START ------------------------------------------------ //"
				break
			end
		a=a+1
		end

		bb=1
		while bb<length(lines)
			if lines[bb]=="#RULES::END -------------------------------------------------- //"
				break
			end
		bb=bb+1
		end

		#Grabbing rules section
		rule_data=lines[a+1:bb-1]

		#Building the reduced rule array for translation to a Dict
		for i in 1:length(rule_data)
			if split(rule_data[i],"=")[2]=="[]"
				rules=push!(rules,"")
			else
				rules=push!(rules,split(rule_data[i],"=")[2])
			end
		end

		#Organizing the information
		Org=Dict("S"=>S,"rxnGeneMat"=>rxnGeneMat,"c"=>c,"b"=>b,"rules"=>rules,"rxns"=>rxns_m,"mets"=>mets,"rev"=>rev,"lb"=>lb,"ub"=>ub,"genes"=>genes,"subSystems"=>subsystems,"rxnNames"=>rxn_des,"metNames"=>metNames)
		
		#Compressing data
		Compress=Dict(model_name=>Org)

		#Saving the data as a .mat file
		matwrite(New_file_name,Compress)
	end

end
