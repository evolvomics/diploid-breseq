/*****************************************************************************
 
 AUTHORS
 
 Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 David B. Knoester
 
 LICENSE AND COPYRIGHT
 
 Copyright (c) 2010 Michigan State University
 
 breseq is free software; you can redistribute it and/or modify it under the  
 terms the GNU General Public License as published by the Free Software 
 Foundation; either version 1, or (at your option) any later version.
 
 *****************************************************************************/

#include "libbreseq/annotated_sequence.h"

#include "libbreseq/error_count.h"

using namespace std;

namespace breseq {
  
  const string BULL_DUMMY_SEQ_ID = "__BULL_DUMMY_SEQ_ID__";
  
  // Load a complete collection of files and verify that sufficient information was loaded
  void cReferenceSequences::LoadFiles(const vector<string>& file_names)
  {
    list<string> sorted_unique_file_names(file_names.begin(), file_names.end());

    sorted_unique_file_names.unique();//Removes non-unique file names

    for(list<string>::const_iterator it = sorted_unique_file_names.begin(); it != sorted_unique_file_names.end(); it++)
    {
      this->LoadFile(*it);
    }
    
    this->Verify();
    this->m_initialized = true;
  }
  
  void cReferenceSequences::LoadFile(const string& file_name)
  {
    ifstream in(file_name.c_str());
    ASSERT(in.good(), "Could not open reference file: " +file_name);

    //! Step 1: Determine what file format it is from first line
    string first_line;
    getline(in, first_line);
    in.close();

    vector<string> split_first_line = split_on_whitespace(first_line);
    FileType file_type = UNKNOWN;
    // Fasta?
    if (first_line[0] == '>') {
      file_type = FASTA;
    }
    // GenBank?
    else if (first_line.find("LOCUS") != string::npos) {
      file_type = GENBANK;
    }
    // GFF?
    else if (first_line.find("##gff-version 3") != string::npos) {
      file_type = GFF3;
    }
    // BULL? Lines have three "words"
    else if (split_first_line.size() == 3) {
      file_type = BULL;
    }
    // The file was not identified and will most likely cause Breseq to fail.
    else if(file_type == UNKNOWN) {
      ERROR("Could not determine format of reference file: " + file_name);
    }
    //! Step 2: Load appropriate file
    switch (file_type) {
      case GENBANK:
      {
        ReadGenBank(file_name);
      }break;
        
      case FASTA:
      {
        ReadFASTA(file_name);
      }break;
        
      case BULL:
      {
        ReadBull(file_name);
      }break;
        
      case GFF3:
      {
        ReadGFF(file_name);
      }break;
        
      default:
        WARN("Could not load the reference file: " +file_name);
    }
  }
  
  void cReferenceSequences::Verify()
  {
    bool Error = false;
    stringstream ss;
    for (vector<cAnnotatedSequence>::iterator itr= this->begin(); itr != this->end(); itr++) {
      cAnnotatedSequence& as = *itr;
      if (!as.get_sequence_length()) {
        ss << "No sequence loaded for: " << as.m_seq_id << endl;
        Error = true;
      }
    }
    if (Error) ERROR(ss.str());
    if (this->empty()) ERROR("Reference files were not loaded");
  }

  
  void cReferenceSequences::ReadFeatureTable(const string &file_name) {
    
    ifstream infile(file_name.c_str());
    assert(!infile.fail());
    
    string line;
    getline(infile,line);
    
    while (!infile.eof())
    {      
      if (line[0] == '#') 
      {
        // description line
        if (line.find("##description ") == 0)
        {
          line.erase(0,14);
          size_t pos = line.find(" ");
          if (pos != string::npos) 
          {
            string seq_id = line.substr(0,pos);
            line.erase(0,pos+1);            
            uint32_t seq_idx = seq_id_to_index(seq_id);
            (*this)[seq_idx].m_description = line;
          }
        }
      }
      else 
      {
        // split line on tabs
        char * cstr = new char [line.size()+1];
        strcpy (cstr, line.c_str());
        
        cSequenceFeaturePtr fp(new cSequenceFeature);
        cSequenceFeature& feature = *fp;
        
        string seq_id;
        
        char * pch;
        pch = strtok(cstr,"\t");
        uint32_t on_token = 0;
        while (pch != NULL)
        {
          //printf ("%s\n",pch);
          
          switch (on_token) {
            case 0:
              seq_id = pch;
              break;
              
            case 1:
              feature["type"] = pch;
              break;
              
            case 2:
              feature["accession"] = pch;
              break;
              
            case 3:
              feature["name"] = pch;
              break;
              
            case 4:
              feature.m_start = atoi(pch);
              break;
              
            case 5:
              feature.m_end = atoi(pch);
              break;
              
            case 6:
              feature.m_strand = atoi(pch);
              break;
              
            case 7:
              feature["product"] = pch;
              break;
          }
          
          pch = strtok (NULL, "\t");
          on_token++;
        }  
        
        delete[] cstr;
        
        (*this)[seq_id].add_feature(fp);
      }
      getline(infile,line);
    }
    
  }
  
  void cReferenceSequences::WriteFeatureTable(const std::string &file_name) {

    ofstream out(file_name.c_str());
    assert(!out.fail()); 
    
    out << "##gff-version\t3" << endl;
    for(vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
    
      
      uint32_t m_length;
      string m_definition, m_version, m_seq_id;
      
      out << "##sequence-region " << it_as->m_seq_id << " 1 " << it_as->m_length << endl;
      out << "##description " << it_as->m_seq_id << " " << it_as->m_description << endl;

    }
    
    for(vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
      
      for (vector<cSequenceFeaturePtr>::iterator it = it_as->m_features.begin(); it < it_as->m_features.end(); it++) {
        cSequenceFeature& feat = **it;
        
        out << it_as->m_seq_id;
        
        if (feat.SafeGet("type") == "") 
          out << "\t" << ".";
        else
          out << "\t" << feat["type"];
        
        if( feat.SafeGet("accession") == "" )
          out << "\t" << ".";
        else
          out << "\t" <<feat["accession"];
        
        if( feat.SafeGet("name") == "" )
          out << "\t" << ".";
        else
          out << "\t" <<feat["name"];
        
        out << "\t" << feat.m_start;
        out << "\t" << feat.m_end;
        out << "\t" << (int)feat.m_strand;
        
        if( feat.SafeGet("product") == "" )
          out << "\t" << ".";
        else
          out << "\t" <<feat["product"];
        out << endl;
      }
      
    }
  }
  
  void cReferenceSequences::ReadFASTA(const std::string &file_name) {
    
    cFastaFile ff(file_name, ios_base::in);
    cFastaSequence on_seq;

    uint32_t on_seq_id = this->size();

    while ( ff.read_sequence(on_seq) ) {
      
      // @JEB sorting the input files by type seems like a better way of doing this,
      //    but we can't do it just by looking at file name endings...
      //
      // If we have a dummy sequence id from loading a BULL feature file before a FASTA,
      // then fill it in and update the index
      if (m_seq_id_to_index.count(BULL_DUMMY_SEQ_ID))
      {
        uint32_t seq_index = this->m_seq_id_to_index[BULL_DUMMY_SEQ_ID];
        this->m_seq_id_to_index.erase(BULL_DUMMY_SEQ_ID);
        (*this)[seq_index].m_seq_id = on_seq.m_name;
        m_seq_id_to_index[(*this)[seq_index].m_seq_id] = seq_index;       
      }
      else
      {
        this->add_new_seq(on_seq.m_name);
      }
      
      // copy the info over (could define an assignment operator...)
      cAnnotatedSequence& this_seq = (*this)[on_seq.m_name];
      this_seq.m_fasta_sequence = on_seq;
      this_seq.m_seq_id = on_seq.m_name;
      this_seq.m_description = on_seq.m_description;
      this_seq.m_length = on_seq.m_sequence.size();
    }


  }
    
  void cReferenceSequences::ReadFASTA(cFastaFile& ff) {
    cFastaSequence on_seq;

    while ( ff.read_sequence(on_seq) ) {
      this->add_new_seq(on_seq.m_name);
      cAnnotatedSequence& new_seq = (*this)[on_seq.m_name];
      new_seq.m_fasta_sequence = on_seq;
      new_seq.m_seq_id = on_seq.m_name;
      new_seq.m_length = on_seq.m_sequence.size();
    }

  }

  void cReferenceSequences::WriteFASTA(const std::string &file_name) {
    
    cFastaFile ff(file_name, ios_base::out);
    for(vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
      ff.write_sequence(it_as->m_fasta_sequence);
    }
  }

  void cReferenceSequences::WriteFASTA(cFastaFile& ff) {

    for(vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
      ff.write_sequence(it_as->m_fasta_sequence);
    }
  }

  /*! ReadGFF abides by the following format:
    http://www.sequenceontology.org/gff3.shtml !*/
  void cReferenceSequences::ReadGFF(const string& file_name)
  {
    cFastaFile in(file_name.c_str(), ios_base::in);
    ASSERT(!in.fail(), "Could not open GFF file: " + file_name);
    //! Step 1: Header //! cFastaFile removes first line, so no header...
      string line;
      while (!in.eof() && getline(in,line)) {
        cSequenceFeaturePtr fp(new cSequenceFeature);
        cSequenceFeature& feature = *fp;

    //! Step 2: Check for GFF Tags, reads FASTA from here.
        // We are concerned about a couple of ## GFF tags
        if (line[0] == '#') {
          if (line.find("##species") != string::npos) {
            //TODO @GRC
            continue;
          }
          else if (line.find("##sequence-region") != string::npos) {
            // Line of form ##sequence-region seqid start end
            stringstream ls(line);
            string x, seq_id, start, end;
            ls >> x >> seq_id >> start >> end;
            this->add_new_seq(seq_id);
            (*this)[seq_id].m_length = from_string<uint32_t>(end);
            continue;
          }
          // Find embedded fasta file
          else if (line.find("##FASTA") != string::npos) {
            /* Things admittedly get a bit hairy right here, you have to
              take the next line (The one with ">XXXXXX") and set it as the
              current line for cFastaFile, after you've read the fasta you
              can exit the function since nothing should exist after.
              */
            getline(in,line);
            in.set_current_line(line);
            this->ReadFASTA(in);
            return; //! Done!
          } else {
            continue;
          }
        }
    /*! Step 3: Split line on tabs("\t") until last column, grab last column until endl("\n"),
        the default for getline(). !*/
        
        stringstream ss(line);
        string seq_id, start, end, strand;

        // Handle columns up to the last one, "attributes"
        // Column 1: "seqid"
        getline(ss, seq_id, '\t');
        // Column 2: "source"
        getline(ss, feature["source"], '\t');
        // Column 3: "type"
        getline(ss, feature["type"], '\t');
        // Column 4: "start"
        getline(ss, start, '\t');
        feature.m_start = from_string<uint32_t>(start);
        // Column 5: "end"
        getline(ss, end, '\t');
        feature.m_end = from_string<uint32_t>(end);
        // Column 6: "score"
        getline(ss, feature["score"], '\t');
        // Column 7: "strand"
        getline(ss, strand, '\t');
        feature.m_strand = 0;
        if (strand == "+")
          feature.m_strand = 1;
        else if (strand == "-")
          feature.m_strand = -1;        
        // Column 8: "phase"
        getline(ss, feature["phase"], '\t');
        // Column 9: "attributes"
        string raw_attributes;

        // Handle parsing the attributes
        getline(ss, raw_attributes);
        vector<string> attributes = split(raw_attributes, ";");

        //Split attribute's key and value by "="
        for (vector<string>::iterator itr = attributes.begin(); itr != attributes.end(); itr++) {
          string& attribute = *itr;
          vector<string> key_value = split(attribute,"=");
          string& key = key_value.front();
          string& value = key_value.back();
        //! Case 2: Multiple values for given key, split by ","
          feature.m_gff_attributes[key] = split(value, ",");
          // unescape special characters after splitting
          for (uint32_t i=0; i<feature.m_gff_attributes[key].size(); i++)
            feature.m_gff_attributes[key][i] = GFF3UnescapeString(feature.m_gff_attributes[key][i]);
        }
              
        // Load certain information into the main hash, so breseq knows to use it
        if (feature.m_gff_attributes.count("Note"))
        {
          feature["product"] = join(feature.m_gff_attributes["Note"], ",");
        }
        
        if (feature.m_gff_attributes.count("Alias"))
        {
          feature["accession"] = join(feature.m_gff_attributes["Alias"], ",");
        }

        if (feature.m_gff_attributes.count("Name"))
        {
          feature["name"] = join(feature.m_gff_attributes["Name"], ",");
        }
        
      
  //! Step 4: Determine if sequence already exists (find or create if not found)
        this->add_new_seq(seq_id);
        (*this)[seq_id].add_feature(fp);
      
        // If this is a landmark "region" corresponding to the entire chromosome grab extra information
        if ((feature["type"] == "region") && (feature.m_start == 1) && (feature.m_end == (*this)[seq_id].m_length))
        {
          if (feature.m_gff_attributes.count("Is_circular"))
            (*this)[seq_id].m_is_circular = (feature.m_gff_attributes["Is_circular"][0] == "true");
          if (feature.m_gff_attributes.count("Note"))
            (*this)[seq_id].m_description = feature.m_gff_attributes["Note"][0];
        }
      }
  }
  
/*! WriteGFF abides by the following format:
  http://www.sequenceontology.org/gff3.shtml !*/
void cReferenceSequences::WriteGFF( const string &file_name ){

  cFastaFile out(file_name.c_str(), ios_base::out);
  assert(!out.fail());
  //! Step 1: Header
  out << "##gff-version 3" << endl;

  for (vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
    out << "##sequence-region" << "\t" << it_as->m_seq_id << "\t" << "1" << "\t" << it_as->m_length << endl;
  }
  
  //! Step 2: Features
  for (vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
      for (vector<cSequenceFeaturePtr>::iterator it = it_as->m_features.begin(); it != it_as->m_features.end(); it++) {
        cSequenceFeature& feat = **it;

        out << it_as->m_seq_id;

        if (feat.SafeGet("source") == "")
          out << "\t.";
        else
          out << "\t" << feat["source"];
        
        if (feat.SafeGet("type") == "")
          out << "\t.";
        else
          out << "\t" << feat["type"];

        out << "\t" << feat.m_start;
        out << "\t" << feat.m_end;

        if (feat.SafeGet("score") == "")
          out << "\t.";
        else
          out << "\t" << feat["score"];

        if (feat.m_strand > 0)
          out << "\t" << "+";
        else
          out << "\t" << "-";

        if (feat.SafeGet("phase") == "")
          out << "\t" << ".";
        else
          out << "\t" << feat["phase"];

        //Attributes
        out << "\t";
        vector<string> attributes;
        map<string,vector<string> >::const_iterator itr;
        for (itr = feat.m_gff_attributes.begin(); itr != feat.m_gff_attributes.end(); itr++) {
          const string& key = itr->first;
          const vector<string>& values = itr->second;
          
          vector<string> s;
          for (vector<string>::const_iterator it = values.begin(); it != values.end(); it++) {
            s.push_back(GFF3EscapeString(*it));
          }
          
          attributes.push_back(key + "=" + join(s, ","));
        }
        out << join(attributes, ";");

        out << std::endl;
      }
  }
  //! Step 3: Fasta
  out << "##FASTA" << endl;
  for (vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
    this->WriteFASTA(out);
  }
  out.close();
}

void cReferenceSequences::ReadGenBank(const string& in_file_name) {

  ifstream in(in_file_name.c_str(), ios_base::in);
  ASSERT(!in.fail(), "Could not open GenBank file: " + in_file_name);

  uint32_t on_seq_index = this->size();

  while (ReadGenBankFileHeader(in)) {
    
    cAnnotatedSequence& s = this->back();
    
    // add a 'region' feature for GFF3 output
    cSequenceFeaturePtr f(new cSequenceFeature);
    (*f)["type"] = "region";
    f->m_strand = 1;
    f->m_start = 1;
    f->m_end = s.m_length;
    
    if (s.m_is_circular) 
      f->m_gff_attributes["Is_circular"].push_back("true");
    else
      f->m_gff_attributes["Is_circular"].push_back("false");
    
    f->m_gff_attributes["Note"].push_back(s.m_description);
    s.add_feature(f);
    
    ReadGenBankFileSequenceFeatures(in, s);
    ReadGenBankFileSequence(in, s);
    this->set_seq_id_to_index(s.m_seq_id, on_seq_index++);
    
  }
}

bool cReferenceSequences::ReadGenBankFileHeader(ifstream& in) {

  //std::cout << "header" << std::endl;
  string line;
  bool found_LOCUS_line = false;
  cAnnotatedSequence* s = NULL;
  while (!in.eof()) {
    std::getline(in, line);
    string first_word = GetWord(line);
    RemoveLeadingTrailingWhitespace(line);

    // This is the first line
    if (first_word == "LOCUS") {

      // Example line 
      // LOCUS       EB03                 4629813 bp    DNA     circular BCT 13-APR-2007
      
      string w;
      w = GetWord(line);
      string seq_id = w;
      
      this->add_new_seq(seq_id);
      s = &((*this)[seq_id]);
      
      w = GetWord(line);
      s->m_length = atoi(w.c_str());
      
      w = GetWord(line);
      w = GetWord(line);
      w = GetWord(line);
      if (to_lower(w) == "circular")
        s->m_is_circular = true;

      // Should only be one line like this per record!
      ASSERT(!found_LOCUS_line, "Multiple LOCUS lines found in single GenBank record."); 
      found_LOCUS_line = true;
    }

    if (first_word == "DEFINITION") {
      ASSERT(s, "Missing LOCUS line before DEFINITION line in GenBank record.");
      s->m_description = line;
    }

    if (first_word == "FEATURES") break;
  }

  return (found_LOCUS_line);
}


void cSequenceFeature::ReadGenBankCoords(string& s, ifstream& in) {

  //std::cerr << "whole: " << s << std::endl;

  // Read through all parentheses
  int32_t parentheses_level = 0;
  size_t parenthesis_pos = s.find_first_of("()");

  while(parenthesis_pos != string::npos) {
    if (s.at(parenthesis_pos) == '(') {
      parentheses_level++;
    } else {
      parentheses_level--;
    }
    parenthesis_pos = s.find_first_of("()", parenthesis_pos+1);
  }

  string value = s;

  // Multiline
  while ((parentheses_level > 0) && !in.eof()) {
    std::getline(in, s);
    RemoveLeadingTrailingWhitespace(s);

    size_t on_pos = s.find_first_of("()");
    while(on_pos != string::npos) {

      if (s.at(on_pos) == '(') {
        parentheses_level++;
      } else {
        parentheses_level--;
      }
      on_pos = s.find_first_of("()", on_pos+1);
    }
    value += s;
  }
  assert(parentheses_level == 0);

  m_strand = 1;
  size_t start_complement = value.find("complement(");

  if (start_complement != string::npos) {
    m_strand = -1;
  }

  size_t start_start_pos = value.find_first_of("1234567890");
  assert(start_start_pos != string::npos);
  size_t start_end_pos = value.find_first_not_of("1234567890", start_start_pos);
  if (start_end_pos == string::npos) start_end_pos = value.size() ;
  string start_s = value.substr(start_start_pos,start_end_pos-start_start_pos);

  size_t end_end_pos = value.find_last_of("1234567890");
  assert(end_end_pos != string::npos);
  size_t end_start_pos = value.find_last_not_of("1234567890", end_end_pos);
  if (end_start_pos == string::npos) start_end_pos = -1;

  string end_s = value.substr(end_start_pos+1,end_end_pos+1 - (end_start_pos+1));

  m_start = atoi(start_s.c_str());
  m_end = atoi(end_s.c_str());
}

void cSequenceFeature::ReadGenBankTag(std::string& tag, std::string& s, std::ifstream& in) {

  // delete leading slash
  tag.erase(0,1);

  // erase through the equals on this line
  int pos = s.find("=");
  s.erase(0,pos+1);

  // If there is a quote then we need to continue to the last quote
  size_t first_quote_pos = s.find("\"");

  if (first_quote_pos == string::npos) {
    (*this)[tag] = s;
    return;
  }

  s.erase(0,first_quote_pos+1);

  size_t second_quote_pos = s.find("\"");

  // One liner
  if (second_quote_pos != string::npos) {
    s.erase(second_quote_pos,s.length());
    (*this)[tag] = s;
    return;
  }

  // If the value is quoted, we have to read additional lines until end quote
  string value = s;

  bool found_last_quote = false;
  while (!found_last_quote && !in.eof()) {
    std::getline(in, s);
    RemoveLeadingTrailingWhitespace(s);

    second_quote_pos = s.find("\"");
    if (second_quote_pos != string::npos) {
      s.erase(second_quote_pos,s.length());
      found_last_quote = true;
    }

    if (tag != "translation") value += " ";
    value += s;
  }
  assert(found_last_quote);

  (*this)[tag] = value;
}

void cReferenceSequences::ReadGenBankFileSequenceFeatures(std::ifstream& in, cAnnotatedSequence& s) {
  //std::cout << "features" << std::endl;
  cSequenceFeature* current_feature(NULL);
  cSequenceFeatureList all_features; // make preliminary list then add once entries are complete
  string line;
  while (!in.eof()) {
    getline(in, line);
    //debug
    //cout << line << endl;
    string first_word = GetWord(line);

    //std::cout << first_word << "::" << line << std::endl;

    // Done with this section...
    if (first_word == "ORIGIN") break;

    // Major tag = new feature or information block
    if (first_word[0] != '/') {

      if (first_word != "BASE") {
        cSequenceFeaturePtr new_feat(new cSequenceFeature);
        all_features.push_back( new_feat );
        current_feature = new_feat.get();
        (*current_feature)["type"] = first_word;
        // parse the rest of the line
        std::string coord_s = GetWord(line);
        current_feature->ReadGenBankCoords(coord_s, in);
      }
    }
    // Minor tag = information about current feature
    else {
      assert(current_feature);

      // Remove leading slash
      current_feature->ReadGenBankTag(first_word, line, in); // reads multi-line entries
    }
  }


  for (cSequenceFeatureList::iterator it = all_features.begin(); it < all_features.end(); it++) {
    cSequenceFeature& feat = **it;

    // common changes for any type
    // use /note as the product if there is not product
    if (feat.SafeGet("product") == "")
    {
      feat["product"] = feat.SafeGet("note");
    }

    if (feat["type"] == "repeat_region") {

      // Don't add unnamed ones to the list...
      //if (it->SafeGet("mobile_element") == "") continue;

      feat["name"] = "repeat_region";

      // E. coli case:
      if (feat.SafeGet("mobile_element") != "")
      {
        feat["name"] = feat["mobile_element"];
        string& name = feat["name"];

        // remove prefix
        int pos = name.find("insertion sequence:");
        if (pos != -1) {
          name.erase(pos,pos+19);
        }

        // remove suffix if "IS(\d)}
        pos = name.find("IS");
        if (pos != -1) {
          int found = name.find_first_not_of("0123456789", 3);
          if (found != -1) {
            name.erase(found, name.length());
          }
        }
      }

      // S. cerevisiae case
      if (feat.SafeGet("rpt_family") != "")
      {
        feat["name"] = feat["rpt_family"];
      }
      //std::cerr << (*it).SafeGet("mobile_element") << " " << (*it).SafeGet("name") << std::endl;

      if (feat.SafeGet("product") == "")
      {
        feat["product"] = "repeat region";
      }
    }
    else
    {
      // Add information
      if (feat.SafeGet("gene") != "") {
        feat["name"] = feat["gene"];
      }
      if ( (feat.SafeGet("name") == "") && (feat.SafeGet("locus_tag") != "") ) {
        feat["name"] = feat["locus_tag"];
      }

      //std::cerr << (*it).SafeGet("name") << " " << (*it).SafeGet("gene") << " " << (*it).SafeGet("locus_tag") << std::endl;

      feat["accession"] = feat.SafeGet("locus_tag");

      // /pseudo tag doesn't take a value
      if (feat.count("pseudo") != 0) {
        feat["type"] = "pseudogene";
      }
    }
    
    s.add_feature(*it);
    
    // transfer to GFF
    if (feat.SafeGet("locus_tag") != "")
      feat.m_gff_attributes["ID"] = make_list<string>(feat["locus_tag"]);
    if (feat.SafeGet("product") != "")
      feat.m_gff_attributes["Note"] = make_list<string>(feat["product"]);
    if (feat.SafeGet("accession") != "")
      feat.m_gff_attributes["Alias"] = make_list<string>(feat["accession"]);
    if (feat.SafeGet("name") != "")
      feat.m_gff_attributes["Name"] = make_list<string>(feat["name"]);
    
    feat["phase"] = "0";
  }
}

void cReferenceSequences::ReadGenBankFileSequence(std::ifstream& in, cAnnotatedSequence& s) {
  //std::cout << "sequence" << std::endl;

  s.m_fasta_sequence.m_name = s.m_seq_id;

  string line;
  while (!in.eof()) {
    getline(in, line);
    string first_word = GetWord(line);
    RemoveLeadingTrailingWhitespace(line);

    //cerr << first_word << "::" << line << std::endl;

    if (first_word == "//") break;
    for(string::iterator i=line.begin(); i!=line.end(); ++i) {
      if (*i == ' ') continue;
      *i = toupper(*i);

      // Scrub nonstandard characters
      if ((*i != 'A' ) && (*i != 'T' ) && (*i != 'C' ) && (*i != 'G' )) {
        *i = 'N';
      }

      s.m_fasta_sequence.m_sequence += *i;
    }
  }

  //cout << s.m_sequence << std::endl;
}

void cReferenceSequences::ReadBull(const string& file_name) {

  // This should only work if we have zero or one reference sequence, because it does not
  // give a seq_id, we assume the FASTA before or the FASTA after is paired with it.
  ASSERT(this->size() <= 1, "Bull format currently only works with one reference sequence");

  // if we encounter the BULL file first, create an empty sequence
  if (this->size() == 0)
  {
    // if one already exists, throw an error
    ASSERT(!this->m_seq_id_to_index.count(BULL_DUMMY_SEQ_ID), 
           "Two BULL formatted feature files loaded in a row. Ambiguous assignment to FASTA files.");
    this->add_new_seq(BULL_DUMMY_SEQ_ID);
  }
    
  ifstream in(file_name.c_str());
  if(!in.good()) WARN("Could not open file:" + file_name);

  char line[10];

  cSequenceFeatureList all_features;

  while ( !in.eof() ) {

    string line;
    getline(in, line);
    vector<string> s = split_on_whitespace(line);

    if (s.size() == 0) continue;
    ASSERT(s.size() == 3, "Expected 3 words in line, found " + to_string(s.size()) + " in line:\n" + line);
    
    cSequenceFeaturePtr current_feature(new cSequenceFeature);
    
    current_feature->m_start = from_string<uint32_t>(s[0]);
    current_feature->m_end = from_string<uint32_t>(s[1]);
    current_feature->m_strand = (current_feature->m_start <= current_feature->m_end) ? 1 : -1;
    
    // transfer to GenBank
    (*current_feature)["name"] = s[2];
    (*current_feature)["type"] = "CDS";
    (*current_feature)["product"] = "";
    (*current_feature)["accession"] = s[2];
    
    // transfer to GFF3
    (*current_feature)["ID"] = s[2];
    (*current_feature)["Note"] = "";
    (*current_feature)["Alias"] = s[2];
    (*current_feature)["Name"] = s[2];
    
    all_features.push_back(current_feature);
  }
  in.close();
  
  cAnnotatedSequence& s = this->back();
  for (cSequenceFeatureList::iterator it = all_features.begin(); it < all_features.end(); it++) {
    cSequenceFeature& feat = **it;
    s.m_features.push_back(*it);
  }
}

string cReferenceSequences::repeat_example(const string &repeat_name, int8_t strand)
{
  vector<cAnnotatedSequence>::iterator itr_seq;

  // loop through all reference sequences
  for (itr_seq = this->begin(); itr_seq != this->end(); itr_seq++) {
    cAnnotatedSequence& this_seq = *itr_seq;
    cSequenceFeatureList& repeats = this_seq.m_repeats;

    cSequenceFeatureList::iterator itr_rep;
    for (itr_rep = repeats.begin(); itr_rep != repeats.end(); itr_rep++) {
      cSequenceFeature& rep = **itr_rep;

      if (rep.SafeGet("name") == repeat_name) {
        string repeat_seq = this_seq.get_sequence_1(rep.m_start, rep.m_end);
        if (strand != rep.m_strand)
          repeat_seq = reverse_complement(repeat_seq);
        return repeat_seq;
      }
    }
  }

  ASSERT(false, "Unknown repeat type: " + repeat_name);
  return "";
}

cSequenceFeature* cReferenceSequences::find_closest_repeat_region_boundary(uint32_t position, vector<cSequenceFeaturePtr>& repeat_list_ref, int32_t max_distance, int32_t direction)
{
  if (repeat_list_ref.size() == 0) return NULL;

  cSequenceFeature* is(NULL);
  int32_t best_distance = max_distance + 1; // this enforces the max distance

  for (uint32_t i = 0; i < repeat_list_ref.size(); i++) //IS
  {
    cSequenceFeature* test_is = repeat_list_ref[i].get();
    
    // Distance from the appropriate end of the repeat
    int32_t test_distance = abs(static_cast<int32_t>(((direction == -1) ? position - test_is->m_end : test_is->m_start - position)));
    
    // We want the closest one without going over that is within max_distance
    if ( (test_distance >= 0) && (test_distance < best_distance) )
    {
      is = test_is;
      best_distance = test_distance;
    }
  }
  return is;
}

cSequenceFeature* cReferenceSequences::get_overlapping_feature(vector<cSequenceFeaturePtr>& feature_list_ref, uint32_t pos)
{
  cSequenceFeature* feat = NULL;
  for (uint32_t i = 0; i < feature_list_ref.size(); i++)
    if (pos >= feature_list_ref[i]->m_start && pos <= feature_list_ref[i]->m_end)
      feat = feature_list_ref[i].get();
  return feat;
}

void cReferenceSequences::find_nearby_genes(
                                            cSequenceFeatureList& gene_list_ref,
                                            uint32_t pos_1,
                                            uint32_t pos_2,
                                            vector<Gene>& within_genes,
                                            vector<Gene>& between_genes,
                                            vector<Gene>& inside_left_genes,
                                            vector<Gene>& inside_right_genes,
                                            Gene& prev_gene,
                                            Gene& next_gene
                                            )
{
  if (pos_2 == UNDEFINED_UINT32)
    pos_2 = pos_1;

  //#	print "$pos_1, $pos_2\n";

  for (uint32_t i = 0; i < gene_list_ref.size(); i++) //GENE
  {
    cSequenceFeature& test_gene_feat = *gene_list_ref[i];
    Gene test_gene = Gene(test_gene_feat); // up-cast, should be a better way to do this @JEB

    if (test_gene.end < pos_1)
      prev_gene = test_gene;

    if (  (test_gene.start <= pos_1) && (test_gene.end >= pos_1)
       && (test_gene.start <= pos_2) && (test_gene.end >= pos_2) )
    {
      within_genes.push_back(test_gene);
      //#	print "^ $test_gene->{name}\n";
    }
    else if ( (test_gene.start <= pos_1) && (test_gene.end >= pos_1) )
    {
      inside_left_genes.push_back(test_gene);
    }
    else if ( (test_gene.start <= pos_2) && (test_gene.end >= pos_2) )
    {
      inside_right_genes.push_back(test_gene);
    }
    else if ( (test_gene.start >= pos_1) && (test_gene.end <= pos_2) )
    {
      between_genes.push_back(test_gene);
      //#	print ">< $test_gene->{name}\n";
    }
    // We've passed the changes, so it is in the previous intergenic space
    if (test_gene.start > pos_2)
    {
      next_gene = test_gene;
      break;
    }
  }

  //#	print "$prev_gene->{name} || $next_gene->{name}\n";

  //return ($prev_gene, $next_gene, \@within_genes, \@between_genes, \@inside_left_genes, \@inside_right_genes);
}

map<string,char> cReferenceSequences::translation_table_11 = make_map<string,char>
  ("TTT", 'F')
  ("TTC", 'F')
  ("TTA", 'L')
  ("TTG", 'L')

  ("TCT", 'S')
  ("TCC", 'S')
  ("TCA", 'S')
  ("TCG", 'S')

  ("TAT", 'Y')
  ("TAC", 'Y')
  ("TAA", '*')
  ("TAG", '*')

  ("TGT", 'C')
  ("TGC", 'C')
  ("TGA", '*')
  ("TGG", 'W')

  ("CTT", 'L')
  ("CTC", 'L')
  ("CTA", 'L')
  ("CTG", 'L')

  ("CCT", 'P')
  ("CCC", 'P')
  ("CCA", 'P')
  ("CCG", 'P')

  ("CAT", 'H')
  ("CAC", 'H')
  ("CAA", 'Q')
  ("CAG", 'Q')

  ("CGT", 'R')
  ("CGC", 'R')
  ("CGA", 'R')
  ("CGG", 'R')

  ("ATT", 'I')
  ("ATC", 'I')
  ("ATA", 'I')
  ("ATG", 'M')

  ("ACT", 'T')
  ("ACC", 'T')
  ("ACA", 'T')
  ("ACG", 'T')

  ("AAT", 'N')
  ("AAC", 'N')
  ("AAA", 'K')
  ("AAG", 'K')

  ("AGT", 'S')
  ("AGC", 'S')
  ("AGA", 'R')
  ("AGG", 'R')

  ("GTT", 'V')
  ("GTC", 'V')
  ("GTA", 'V')
  ("GTG", 'V')

  ("GCT", 'A')
  ("GCC", 'A')
  ("GCA", 'A')
  ("GCG", 'A')

  ("GAT", 'D')
  ("GAC", 'D')
  ("GAA", 'E')
  ("GAG", 'E')

  ("GGT", 'G')
  ("GGC", 'G')
  ("GGA", 'G')
  ("GGG", 'G')
;

char cReferenceSequences::translate(string seq)
{
  return (translation_table_11.count(seq) == 0) ? '?' : translation_table_11[seq];
}

void cReferenceSequences::annotate_1_mutation(diff_entry& mut, uint32_t start, uint32_t end, bool repeat_override)
{
  // this could be moved to the object
  string intergenic_seperator = "/";

  // initialize everything, even though we don"t always use it
  mut["aa_position"] = "";
  mut["aa_ref_seq"] = "";
  mut["aa_new_seq"] = "";
  mut["codon_position"] = "";
  mut["codon_ref_seq"] = "";
  mut["codon_new_seq"] = "";
  mut["gene_name"] = "";
  mut["gene_position"] = "";
  mut["gene_product"] = "";
  mut["gene_list"] = ""; //#affected genes

  string seq_id = mut["seq_id"];

  //or die "Unknown seq_id in reference sequence info: $seq_id\n";

  vector<cSequenceFeaturePtr>& gene_list_ref = (*this)[seq_id].m_genes;
  vector<cSequenceFeaturePtr>& repeat_list_ref = (*this)[seq_id].m_repeats;

  int32_t size = end - start + 1;

  Gene prev_gene, next_gene;
  vector<Gene> within_genes;
  vector<Gene> between_genes;
  vector<Gene> inside_left_genes;
  vector<Gene> inside_right_genes;

  cSequenceFeature* repeat_region(NULL);
  if (repeat_override)
  {
    assert(start == end);
    repeat_region = get_overlapping_feature(repeat_list_ref, start);
    if (repeat_region != NULL)
    {
      Gene within_gene(*repeat_region);
      within_genes.push_back(within_gene);
    }
  }

  if (repeat_region == NULL)
    find_nearby_genes(gene_list_ref, start, end, within_genes, between_genes, inside_left_genes, inside_right_genes, prev_gene, next_gene);

  // Mutation is intergenic
  if (within_genes.size() + between_genes.size() + inside_left_genes.size() + inside_right_genes.size() == 0)
  {
    mut["snp_type"] = "intergenic";

    mut["gene_name"] += (prev_gene.name.size() > 0) ? prev_gene.name : "–"; //en-dash
    mut["gene_name"] += intergenic_seperator;
    mut["gene_name"] += (next_gene.name.size() > 0) ? next_gene.name : "–"; //en-dash

    if (prev_gene.name.size() > 0)
    {
      mut["gene_position"] += "intergenic (";
      mut["gene_position"] += (prev_gene.strand) ? "+" : "-"; //hyphen
      mut["gene_position"] += to_string(start - prev_gene.end);
    }
    else
    {
      mut["gene_position"] += "intergenic (–";
    }
    mut["gene_position"] += intergenic_seperator;
    if (next_gene.name.size() > 0)
    {
      mut["gene_position"] += (next_gene.strand) ? "-" : "+"; //hyphen
      mut["gene_position"] += to_string(next_gene.start - end);
    }
    else
    {
      mut["gene_position"] += "–"; //en-dash
    }
    mut["gene_position"] += ")";

    mut["gene_product"] += (prev_gene.name.size() > 0) ? prev_gene.product : "–"; //en-dash
    mut["gene_product"] += intergenic_seperator;
    mut["gene_product"] += (next_gene.name.size() > 0) ? next_gene.product : "–"; //en-dash

    return;
  }
  // Mutation is completely within genes
  else if (within_genes.size() > 0)
  {
    /// TODO: It can be within multiple genes, in which case we need to annotate
    /// the change it causes in each reading frame UGH! YUCKY!
    /// FOR NOW: just take the first of the within genes...
    Gene gene = within_genes[0];
    mut["gene_name"] = gene.name;
    mut["gene_product"] = gene.product;

    //#added for gene table
    mut["gene_list"] = gene.name;

    int32_t within_gene_start = (gene.strand) ? gene.start : gene.end;

    if (start == end)
    {
      mut["gene_position"] = to_string(abs(static_cast<int32_t>(start) - within_gene_start) + 1);
    }
    else
    {
      uint32_t gene_start = abs(static_cast<int32_t>(start) - within_gene_start) + 1;
      uint32_t gene_end = abs(static_cast<int32_t>(end) - within_gene_start) + 1;
      mut["gene_position"] = (gene_start < gene_end)
        ? to_string(gene_start) + "-" + to_string(gene_end)  //hyphen
        : to_string(gene_end) + "-" + to_string(gene_start); //hyphen
    }

    string gene_nt_size = to_string(gene.end - gene.start + 1);

    // ...but the gene is a pseudogene or not a protein coding gene
    if (gene.pseudogene)
    {
      mut["snp_type"] = "pseudogene";
      mut["gene_position"] = "pseudogene (" + mut["gene_position"] + "/" + gene_nt_size + " nt)";
      return;
    }
    else if (gene.type != "CDS")
    {
      mut["snp_type"] = "noncoding";
      mut["gene_position"] = "noncoding (" + mut["gene_position"] + "/" + gene_nt_size + " nt)";
      return;
    }

    //#only add gene information to SNPs and RA mutations that don"t include indels...
    if ((mut._type != SNP) && !((mut._type == RA) && (mut["ref_base"] != ".") && (mut["new_base"] != ".")))
    {
      mut["gene_position"] = "coding (" + mut["gene_position"] + "/" + gene_nt_size + " nt)";
      return;
    }

    // this is for RA...
    if (!mut.entry_exists("ref_seq")) mut["ref_seq"] = mut["ref_base"];
    if (!mut.entry_exists("new_seq")) mut["new_seq"] = mut["new_base"];

    // determine the old and new translation of this codon
    mut["aa_position"] = to_string((from_string<uint32_t>(mut["gene_position"]) - 1) / 3 + 1); // 1 indexed
    mut["codon_position"] = to_string(int(abs(static_cast<int32_t>(start) - within_gene_start)) % 3 + 1); // 1 indexed

    string& ref_string = (*this)[seq_id].m_fasta_sequence.m_sequence;
    string codon_seq = (gene.strand)
      ? ref_string.substr(gene.start + 3 * (from_string<uint32_t>(mut["aa_position"]) - 1) - 1, 3)
      : reverse_complement(ref_string.substr(gene.end - 3 * from_string<uint32_t>(mut["aa_position"]), 3));

    mut["codon_ref_seq"] = codon_seq;
    mut["aa_ref_seq"] = translate(mut["codon_ref_seq"]);

    mut["codon_new_seq"] = codon_seq;
    //#remember to revcom the change if gene is on opposite strand
    mut["codon_new_seq"][from_string<uint32_t>(mut["codon_position"]) - 1] = gene.strand ?
      mut["new_seq"][0]
      : reverse_complement(mut["new_seq"])[0];
    mut["aa_new_seq"] =  translate(mut["codon_new_seq"]);

    mut["snp_type"] = (mut["aa_ref_seq"] != mut["aa_new_seq"]) ? "nonsynonymous" : "synonymous";
  }

  //The mutation actually contains several genes
  else if (between_genes.size() + inside_left_genes.size() + inside_right_genes.size() > 0)
  {
    vector<Gene> gene_list;
    gene_list.insert( gene_list.end(), inside_left_genes.begin(), inside_left_genes.end() );
    gene_list.insert( gene_list.end(), between_genes.begin(), between_genes.end() );
    gene_list.insert( gene_list.end(), inside_right_genes.begin(), inside_right_genes.end() );

    vector<string> gene_name_list;
    for (vector<Gene>::iterator it=inside_left_genes.begin(); it != inside_left_genes.end(); it++)
    {
      gene_name_list.push_back("<i>[" + it->name + "]</i>");
    }
    for (vector<Gene>::iterator it=between_genes.begin(); it != between_genes.end(); it++)
    {
      gene_name_list.push_back("<i>" + it->name + "</i>");
    }
    for (vector<Gene>::iterator it=inside_right_genes.begin(); it != inside_right_genes.end(); it++)
    {
      gene_name_list.push_back("<i>[" + it->name + "]</i>");
    }

    mut["gene_product"] = join (gene_name_list, ", ");

    if (gene_name_list.size() == 1)
      mut["gene_name"] = gene_name_list[0];
    else
      mut["gene_name"] = gene_name_list.front() + "–" + gene_name_list.back();  //en-dash
  }
}

void cReferenceSequences::annotate_mutations(genome_diff& gd, bool only_muts)
{
  //keep track of other mutations that affect SNPs
  //because we may double-hit a codon

  //TODO: the proper way to do this is to create list of SNPs that have been hit
  // hashed by gene protein accession ID and AA position within gene
  // and have the annotation point to them (and back at them)
  // so that the codon will be correctly updated with all changes and we can notify the
  // changes that their SNP_type is not really SNP, but multiple hit SNP.

  diff_entry_list muts = gd.show_list();
  for (diff_entry_list::iterator it=muts.begin(); it!=muts.end(); it++)
  {
    diff_entry& mut= **it;

    if (only_muts && (mut.is_evidence())) continue;

    if (mut._type == SNP)
    {
      mut["_ref_seq"] = get_sequence_1(mut["seq_id"], from_string<uint32_t>(mut["position"]), from_string<uint32_t>(mut["position"]));
      annotate_1_mutation(mut, from_string<uint32_t>(mut["position"]), from_string<uint32_t>(mut["position"]));
    }
    else if (mut._type == SUB)
    {
      annotate_1_mutation(mut, from_string<uint32_t>(mut["position"]), from_string<uint32_t>(mut["position"]) + from_string<uint32_t>(mut["size"]) - 1);
    }
    else if (mut._type == DEL)
    {
      annotate_1_mutation(mut, from_string<uint32_t>(mut["position"]), from_string<uint32_t>(mut["position"]) + from_string<uint32_t>(mut["size"]) - 1);
    }
    else if (mut._type == INS)
    {
      annotate_1_mutation(mut, from_string<uint32_t>(mut["position"]), from_string<uint32_t>(mut["position"]));
    }
    else if (mut._type == CON)
    {
      annotate_1_mutation(mut, from_string<uint32_t>(mut["position"]), from_string<uint32_t>(mut["position"]) + from_string<uint32_t>(mut["size"]) - 1);
    }
    else if (mut._type == MOB)
    {
      annotate_1_mutation(mut, from_string<uint32_t>(mut["position"]), from_string<uint32_t>(mut["position"]) + from_string<uint32_t>(mut["duplication_size"]) - 1);
    }
    else if (mut._type == INV)
    {
      annotate_1_mutation(mut, from_string<uint32_t>(mut["position"]), from_string<uint32_t>(mut["position"]));
      mut["gene_name_1"] = mut["gene_name"];
      mut["gene_product_1"] = mut["gene_product"];
      annotate_1_mutation(mut, from_string<uint32_t>(mut["position"]) + from_string<uint32_t>(mut["size"])-1, from_string<uint32_t>(mut["position"]) + from_string<uint32_t>(mut["size"])-1);
      mut["gene_name_2"] = mut["gene_name"];
      mut["gene_product_2"] = mut["gene_product"];
      mut._fields.erase("gene_name");
      mut._fields.erase("gene_product");
    }
    else if (mut._type == AMP)
    {
      annotate_1_mutation(mut, from_string<uint32_t>(mut["position"]), from_string<uint32_t>(mut["position"]) + from_string<uint32_t>(mut["size"]) - 1);
    }
    else if (mut._type == JC)
    {
      diff_entry side_1;
      side_1[SEQ_ID] = mut["side_1_seq_id"];
      annotate_1_mutation(side_1, from_string<uint32_t>(mut["side_1_position"]), from_string<uint32_t>(mut["side_1_position"]), true);
      //copy over entries with prefix
      for(diff_entry::map_t::iterator it=side_1._fields.begin(); it!=side_1._fields.end(); it++)
      {
        mut["_side_1"+ it->first] = it->second;
      }

      diff_entry side_2;
      side_2[SEQ_ID] = mut["side_2_seq_id"];
      annotate_1_mutation(side_2, from_string<uint32_t>(mut["side_2_position"]), from_string<uint32_t>(mut["side_2_position"]), true);
      //copy over entries with prefix
      for(diff_entry::map_t::iterator it=side_2._fields.begin(); it!=side_2._fields.end(); it++)
      {
        mut["_side_2"+ it->first] = it->second;
      }
    }
    else if (mut._type == RA)
    {
      annotate_1_mutation(mut, from_string<uint32_t>(mut["position"]), from_string<uint32_t>(mut["position"]));
    }
    else if (mut._type == MC)
    {
      annotate_1_mutation(mut, from_string<uint32_t>(mut["start"]), from_string<uint32_t>(mut["end"]));
    }
  }
}

void cReferenceSequences::polymorphism_statistics(Settings& settings, Summary& summary)
{
  string reference_fasta_file_name = settings.reference_fasta_file_name;
  vector<string> seq_ids = this->seq_ids();

  // some local variable lookups for convenience
  double log10_ref_length = log(this->total_length()) / log(10);

  //
  // Replacement for below
  //
  // ToDo: This should really make a different column for each input read set.
  //
  string coverage_fn = settings.file_name(settings.unique_only_coverage_distribution_file_name, "@", "");
  string outputdir = dirname(coverage_fn) + "/";
  //chomp $outputdir; $outputdir .= "/";
  string count_file_name = outputdir + "error_counts.tab";

  ifstream COUNT(count_file_name.c_str());
  assert(COUNT.is_open());
  string count_header_line, count_covariate_line;
  getline(COUNT, count_covariate_line); // ignore the first line

  // we parse the covariates to know how many different qualities we are dealing with...
  covariates_used_t         covariate_used;
  covariates_max_t          covariate_max;
  covariates_enforce_max_t  covariate_enforce_max;
  covariates_offset_t       covariate_offset;
  bool                      per_position;
  cErrorTable::read_covariates(count_covariate_line, covariate_used, covariate_max, covariate_enforce_max, covariate_offset, per_position);

  getline(COUNT, count_header_line); // second line is the header

  vector<string> count_header_list = split(count_header_line, "\t");

  uint32_t count_column = UNDEFINED_UINT32;
  uint32_t quality_column = UNDEFINED_UINT32;
  for (uint32_t i = 0; i < count_header_list.size(); i++)
  {
    if (count_header_list[i] == "quality")
      quality_column = i;
    else if (count_header_list[i] == "count")
      count_column = i;
  }

  ASSERT( (quality_column != UNDEFINED_UINT32) && (count_column != UNDEFINED_UINT32),
          "'quality' and 'count' columns not found in file: " + count_file_name);


  ////#print "$count_column $quality_column\n";

  vector<uint32_t> quality_count_list(covariate_max[k_quality]);
  string line;
  while (COUNT.good())
  {
    getline(COUNT, line);
    vector<string> line_list = split(line,  "\t");
    if (line_list.size() < 2) break; // empty line
    uint32_t count = from_string<uint32_t>(line_list[count_column]);
    uint32_t quality = from_string<uint32_t>(line_list[quality_column]);
    quality_count_list[quality] += count;
  }
  COUNT.close();

  string genome_error_counts_file_name = settings.genome_error_counts_file_name;

  ofstream GEC(genome_error_counts_file_name.c_str());
  assert(GEC.is_open());
  for (uint32_t i = 1; i < quality_count_list.size(); i++)
  {
    uint32_t val = 0;
    //if (is_defined(quality_count_list[i]))
      val = quality_count_list[i];
    GEC << val << endl;
  }
  GEC.close();

  string polymorphism_statistics_input_file_name = settings.polymorphism_statistics_input_file_name;
  string polymorphism_statistics_output_file_name = settings.polymorphism_statistics_output_file_name;

  /// Load the older GenomeDiff and add new fields
  string ra_mc_genome_diff_file_name = settings.ra_mc_genome_diff_file_name;
  genome_diff gd(ra_mc_genome_diff_file_name);

  string polymorphism_statistics_r_script_file_name = settings.polymorphism_statistics_r_script_file_name;
  string polymorphism_statistics_r_script_log_file_name = settings.polymorphism_statistics_r_script_log_file_name;
  uint32_t total_reference_length = summary.sequence_conversion.total_reference_sequence_length;

  string command = "R --vanilla total_length=" + to_string<uint32_t>(total_reference_length) + " in_file=" + polymorphism_statistics_input_file_name + " out_file=" + polymorphism_statistics_output_file_name + " qual_file=" + genome_error_counts_file_name + " < " + polymorphism_statistics_r_script_file_name + " > " + polymorphism_statistics_r_script_log_file_name;
  SYSTEM(command);

  // Read R file and add new results corresponding to all columns
  ifstream ROUT(polymorphism_statistics_output_file_name.c_str());
  assert(ROUT.is_open()); // or die "Could not find file: $polymorphism_statistics_output_file_name";
  string header;
  getline(ROUT, header);
  //chomp $header;
  vector<string> header_list = split(header, "\t");

  genome_diff new_gd;

  diff_entry_list muts = gd.evidence_list();
  for (diff_entry_list::iterator it=muts.begin(); it!=muts.end(); it++)
  {
    diff_entry& mut= **it;

    // lines only exist for RA evidence
    if (mut._type != RA)
    {
      new_gd.add(mut);
      continue;
    }

    // lines only exist for polymorphisms
    if ((mut[FREQUENCY] == "1") || (mut[FREQUENCY] == "0"))
    {
      new_gd.add(mut);
      continue;
    }

    string line;
    getline(ROUT, line);
    //chomp $header;
    vector<string> line_list = split(line, "\t");

    for (uint32_t j = 0; j < header_list.size(); j++)
    {
      assert(line_list.size() > j); // die "Incorrect number of items on line:\n$line"
      mut[header_list[j]] = line_list[j];
    }

    // Evalue cutoff again (in case we are only running this part)
    if (from_string<double>(mut["polymorphism_quality"]) < settings.polymorphism_log10_e_value_cutoff)
      add_reject_reason(mut, "EVALUE");

    // Frequency cutoff
    if ( (from_string<double>(mut["frequency"]) < settings.polymorphism_frequency_cutoff)
      || (from_string<double>(mut["frequency"]) > 1-settings.polymorphism_frequency_cutoff) )
      add_reject_reason(mut, "POLYMORPHISM_FREQUENCY_CUTOFF");

    // Minimum coverage on both strands
    double polymorphism_coverage_limit_both_bases = settings.polymorphism_coverage_both_strands;
    bool passed = true;
    vector<string> top_bot = split(mut["ref_cov"], "/");
    double top = from_string<double>(top_bot[0]);
    double bot = from_string<double>(top_bot[1]);
    passed = passed && (top >= polymorphism_coverage_limit_both_bases);
    top_bot = split(mut["new_cov"], "/");
    top = from_string<double>(top_bot[0]);
    bot = from_string<double>(top_bot[1]);
    passed = passed && (top >= polymorphism_coverage_limit_both_bases);
    passed = passed && (bot >= polymorphism_coverage_limit_both_bases);

    if (!passed)
      add_reject_reason(mut, "POLYMORPHISM_STRAND");
    if (from_string<double>(mut["ks_quality_p_value"]) < settings.polymorphism_bias_p_value_cutoff)
      add_reject_reason(mut, "KS_QUALITY_P_VALUE");
    if (from_string<double>(mut["fisher_strand_p_value"]) < settings.polymorphism_bias_p_value_cutoff)
      add_reject_reason(mut, "FISHER_STRAND_P_VALUE");

    ////// Optionally, ignore if in a homopolymer stretch
    if (settings.polymorphism_reject_homopolymer_length != UNDEFINED_UINT32)
    {
      uint32_t test_length = 20;
      string seq_id = mut["seq_id"];
      uint32_t end_pos = from_string<uint32_t>(mut["position"]);
      uint32_t start_pos = end_pos - test_length + 1;
      if (start_pos < 1) start_pos = 1;
      string bases = this->get_sequence_1(seq_id, start_pos, end_pos);

      //#print Dumper($mut);
      //#print "$bases\n";

      uint32_t same_base_length = 0;
      string first_base = bases.substr(end_pos - start_pos, 1);
      for (uint32_t j = end_pos; j >= start_pos; j--)
      {
        string this_base = bases.substr(j - start_pos, 1);
        if (first_base != this_base) break;
        same_base_length++;
      }

      //#print "$same_base_length\n";
      if (same_base_length >= settings.polymorphism_reject_homopolymer_length)
      {
        add_reject_reason(mut, "HOMOPOLYMER_STRETCH");
      }
    }

    if (
        mut.number_reject_reasons() > 0
        && (from_string<double>(mut[POLYMORPHISM_QUALITY]) > settings.mutation_log10_e_value_cutoff)
        && (from_string<double>(mut[FREQUENCY]) > 0.5)
        )
    {
      //#print Dumper($mut);
      mut[FREQUENCY] = "1";
      mut._fields.erase(REJECT);

      // FIX -- need to re-evaluate whether it would have been accepted as a normal mutation
      // This is NOT the right quality being used here. Need a separate quality for consensus call and polymorphism call!
      if (from_string<double>(mut[POLYMORPHISM_QUALITY]) < settings.mutation_log10_e_value_cutoff)
        add_reject_reason(mut, "EVALUE");
    }

    new_gd.add(mut);

    // END EXPERIMENTAL
  }

  ROUT.close();

  /// Write out the file which now has much more data
  string polymorphism_statistics_ra_mc_genome_diff_file_name = settings.polymorphism_statistics_ra_mc_genome_diff_file_name;
  new_gd.write(polymorphism_statistics_ra_mc_genome_diff_file_name);

}

std::string GetWord(std::string &s) {
  RemoveLeadingWhitespace(s);
  int found = s.find_first_of(" =");
  std::string w = s.substr(0,found);
  s.erase(0, found);
  return w;
}


void RemoveLeadingWhitespace(std::string &s) {
  int found = s.find_first_not_of(" \t");
  s.erase(0,found);
}

void RemoveLeadingTrailingWhitespace(std::string &s) {
  int found = s.find_first_not_of(" \t\n\r");
  s.erase(0,found);
  found = s.find_last_not_of(" \t\n\r");
  s.erase(found+1,s.length());
}

uint32_t alignment_mismatches(const alignment_wrapper& a, const cReferenceSequences& ref_seq_info)
{
  bool verbose = false;
  uint32_t mismatches = 0;

  uint32_t index = a.reference_target_id();
  const string& const_ref_string = ref_seq_info[index].m_fasta_sequence.m_sequence;
  string ref_string = const_ref_string.substr(a.reference_start_0(), a.reference_match_length());
  uint32_t ref_pos = 0;

  string read_string = a.read_char_sequence().substr(a.query_start_0(), a.query_match_length());
  uint32_t read_pos = 0;

  uint32_t* cigar_list = a.cigar_array(); // cigar array for this alignment

  if (verbose)
  {
    cout << a.read_name() << endl;
  }

  for (uint32_t i = 0; i < a.cigar_array_length(); i++)
  {
    char op = cigar_list[i] & BAM_CIGAR_MASK;
    uint32_t len = cigar_list[i] >> BAM_CIGAR_SHIFT;

    // soft padding counts as a mismatch
    if (op == BAM_CSOFT_CLIP)
    {
      mismatches += len;
    }
    else if (op == BAM_CDEL)
    {
      mismatches += len;
      ref_pos += len;
    }
    else if (op == BAM_CINS)
    {
      mismatches += len;
      read_pos += len;
    }
    else if (op == BAM_CMATCH)
    {
      for (uint32_t j = 0; j < len; j++)
      {
        if (verbose)
        {
          cout << "REF: " << ref_pos << "  " << ref_string[ref_pos] << endl;
          cout << "READ: " << read_pos << "  " << read_string[read_pos] << endl;
        }
        if (ref_string[ref_pos] != read_string[read_pos])
        {
          mismatches++;
        }
        read_pos++;
        ref_pos++;
      }
    }
  }

  //#	print $a->qname . "\n$mismatches\n$cigar_string\n$ref_string\n$read_string\n" if ($mismatches);
  return mismatches;
}

string shifted_cigar_string(const alignment_wrapper& a, const cReferenceSequences& ref_seq_info)
{
  bool verbose = false;

  string ref_seq = ref_seq_info[a.reference_target_id()].m_fasta_sequence.m_sequence;
  uint32_t ref_seq_index = 0;
  string read_seq = a.read_char_sequence();
  uint32_t read_seq_index = 0;
  vector<pair<char,uint16_t> > cigar_pair_array = a.cigar_pair_array();

  if (verbose)
  {
    cout << a.read_name() << endl;

    if (a.read_name() == "GW1ULQG02EIUY7")
    {
      cout << "debug" << endl;
    }
    cout << a.cigar_string() << endl;
  }

  for (vector<pair<char,uint16_t> >::iterator it = cigar_pair_array.begin(); it != cigar_pair_array.end(); it++)
  {
    char op = it->first;
    uint16_t len = it->second;

    if (op == 'D')
    {
      char base = ref_seq[ref_seq_index];
      bool all_same_base = true;
      for (uint32_t j = 1; j < len; j++)
        all_same_base = all_same_base && (ref_seq[ref_seq_index + j] == base);

      if (all_same_base)
      {
        uint32_t shift_amount = 0;
        while (ref_seq[ref_seq_index + len + shift_amount] == base)
          shift_amount++;

        if (shift_amount > 0)
        {
          if (verbose)
          {
            cout << "Shifting D: " << a.read_name() << endl;
          }

          if (it != cigar_pair_array.begin())
            (it - 1)->second += shift_amount;
          else
            cigar_pair_array.insert(it, make_pair<char,uint16_t>('M', shift_amount));

          if (it != cigar_pair_array.end())
            (it + 1)->second -= shift_amount;
          else
            cigar_pair_array.push_back(make_pair<char,uint16_t>('M', shift_amount));
        }

        ref_seq_index += len + shift_amount;
      }
      else
      {
        ref_seq_index += len;
      }
    }
    else if (op == 'I')
    {
      char base = read_seq[read_seq_index];
      bool all_same_base = true;
      for (uint32_t j = 1; j < len; j++)
        all_same_base = all_same_base && (read_seq[read_seq_index + j] == base);

      if (all_same_base)
      {
        uint32_t shift_amount = 0;
        while (read_seq[read_seq_index + len + shift_amount] == base)
          shift_amount++;

        if (shift_amount > 0)
        {
          if (verbose)
          {
            cout << "Shifting I: " << a.read_name() << endl;
          }

          if (it != cigar_pair_array.begin())
            (it - 1)->second += shift_amount;
          else
            cigar_pair_array.insert(it, make_pair<char,uint16_t>('M', shift_amount));

          if (it != cigar_pair_array.end())
            (it + 1)->second -= shift_amount;
          else
            cigar_pair_array.push_back(make_pair<char,uint16_t>('M', shift_amount));
        }

        read_seq_index += len + shift_amount;
      }
      else
      {
        read_seq_index += len;
      }
    }
    else
    {
      ref_seq_index += len;
      read_seq_index += len;
    }
  }

  stringstream shifted_cigar_string_ss;
  for (vector<pair<char,uint16_t> >::iterator it = cigar_pair_array.begin(); it != cigar_pair_array.end(); it++)
    shifted_cigar_string_ss << it->second << it->first;

  string shifted_cigar_string = shifted_cigar_string_ss.str();
  if (verbose)
  {
    cout << shifted_cigar_string << endl;
  }

  return shifted_cigar_string;
}

} // breseq namespace

