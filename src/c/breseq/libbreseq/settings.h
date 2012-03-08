/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

 *****************************************************************************/


#ifndef _BRESEQ_SETTINGS_H_
#define _BRESEQ_SETTINGS_H_

#include "common.h"
#include "storable.h"


namespace breseq
{

	class ExecutionTime : public Storable {
  public:
		string _message;
		string _name;
		time_t _time;
		string _formatted_time;
		double _time_elapsed; // in seconds
		string _formatted_time_elapsed;
		time_t _time_start;
		string _formatted_time_start;
		time_t _time_end;
		string _formatted_time_end;
    
    void serialize(ofstream& f)
    {
      write_to_file(f, _message);
      write_to_file(f, _name);
      write_to_file(f, _time);
      write_to_file(f, _formatted_time);
      write_to_file(f, _time_elapsed);
      write_to_file(f, _formatted_time_elapsed);
      write_to_file(f, _time_start);
      write_to_file(f, _formatted_time_start);
      write_to_file(f, _time_end);
      write_to_file(f, _formatted_time_end);
    }
    void deserialize(ifstream& f)
    {
      read_from_file(f, _message);
      read_from_file(f, _name);
      read_from_file(f, _time);
      read_from_file(f, _formatted_time);
      read_from_file(f, _time_elapsed);
      read_from_file(f, _formatted_time_elapsed);
      read_from_file(f, _time_start);
      read_from_file(f, _formatted_time_start);
      read_from_file(f, _time_end);
      read_from_file(f, _formatted_time_end);
    }
	};

	class Coverage : public Storable
	{
			public:
				double deletion_coverage_propagation_cutoff;
				double deletion_coverage_seed_cutoff;
        double junction_coverage_cutoff;
				double nbinom_size_parameter;
				double nbinom_mean_parameter;
				double nbinom_prob_parameter;
				double average;
				double variance;
				double dispersion;
    
    void serialize(ofstream& f)
    {
      write_to_file(f, deletion_coverage_propagation_cutoff);
      write_to_file(f, deletion_coverage_seed_cutoff);
      write_to_file(f, junction_coverage_cutoff);
      write_to_file(f, nbinom_size_parameter);
      write_to_file(f, nbinom_mean_parameter);
      write_to_file(f, nbinom_prob_parameter);
      write_to_file(f, average);
      write_to_file(f, variance);
      write_to_file(f, dispersion);
    }
    void deserialize(ifstream& f)
    {
      read_from_file(f, deletion_coverage_propagation_cutoff);
      read_from_file(f, deletion_coverage_seed_cutoff);
      read_from_file(f, junction_coverage_cutoff);
      read_from_file(f, nbinom_size_parameter);
      read_from_file(f, nbinom_mean_parameter);
      read_from_file(f, nbinom_prob_parameter);
      read_from_file(f, average);
      read_from_file(f, variance);
      read_from_file(f, dispersion);
    }
	};
	
	class cReferenceSequences;

	// We need to be able to group read files for two reasons
	// 1) They may be paired-end, so we want to map them together
	// 2) They may have the same error rates, so we want to treat them together for error analysis

	struct cReadFile
	{
	public:
		string m_original_file_name;  // the original name provided at the command line
		string m_base_name;           // the original name minus path and .fastq ending (if any)
    string m_converted_file_name; // the name of the converted FASTQ file (if it exists)
		uint32_t m_paired_end_group;  // indicates what file contains paired reads
		uint32_t m_error_group;       // indicates what other read files have the same error rates
		uint32_t m_id;                // index used to refer to this fastq file in BAM
    
    cReadFile()
    {
      m_paired_end_group = UINT_MAX;
      m_error_group = UINT_MAX;
      m_id = UINT_MAX;
    }
    
    string file_name()
    {
      if (m_converted_file_name != "") return m_converted_file_name; 
      return m_original_file_name;
    }
    
    string base_name() { return m_base_name; }
	};

	class cReadFiles : public vector<cReadFile>
	{
	public:
    map<string,string> read_file_to_fastq_file_name_map;
    map<string,string> read_file_to_converted_fastq_file_name_map;

    
		cReadFiles() { };
		cReadFiles(const vector<string>& read_file_names) { Init(read_file_names); };
		~cReadFiles() { };

		void Init(const vector<string>& read_file_names);
    string base_name_to_read_file_name(const string& base_name);
    vector<string> base_names()
    {
      vector<string> return_value;
      for(vector<cReadFile>::iterator it=this->begin(); it!=this->end(); it++)
      {
        return_value.push_back(it->base_name());
      }
      return return_value;
    }
    
	};

	struct Settings
	{
	public:
    
    ////////////////////
    //! Data
    ////////////////////
    
    string byline;
		string website;
    string full_command_line;
		string arguments;
    
    //! Done file tracking
    map<string,string> done_key_messages;
		vector<ExecutionTime> execution_times;
    
    //! Read files
    cReadFiles read_files;
    
    ////////////////////
    //! Settings
    ////////////////////
    
    //! Settings: Global Workflow and Output
    string base_output_path;              // Default = cwd COMMAND-LINE OPTION
    vector<string> read_file_names;       // REQUIRED COMMAND-LINE OPTION
    vector<string> reference_file_names;  // REQUIRED COMMAND-LINE OPTION
    string run_name;          // Default = <none> COMMAND-LINE OPTION
    string print_run_name;    // run_name with '_' replaced by ' '
    
    //! Options that control which parts of the pipeline to execute
    bool no_junction_prediction;      // Default = false COMMAND-LINE OPTION
		bool no_mutation_prediction;      // Default = false
		bool no_deletion_prediction;      // Default = false
		bool no_alignment_generation;     // Default = false
		bool do_copy_number_variation;    // Default = false COMMAND-LINE OPTION
		bool do_periodicity;							// Default = false COMMAND-LINE OPTION

    //! DEBUG options
    
    // verbose level
    uint32_t verbose;                         // Default = 0 (OFF) COMMAND-LINE OPTION
		uint32_t alignment_read_limit;            // Default = 0 (OFF)
		uint32_t candidate_junction_read_limit;   // Default = 0 (OFF)
    uint32_t resolve_alignment_read_limit;    // Default = 0 (OFF)
    //! Output unmatched read file?
		bool no_unmatched_reads;                  // Default = false
    //! Don't delete intermediate files
    bool keep_all_intermediates;              // Default = false

    //! Settings: Read Alignment and Candidate Junction Read Alignment
    uint32_t ssaha2_seed_length;  // Default = 13
    uint32_t ssaha2_skip_length;  // Default = 1 (i.e. no skipping)
    
    //! reads are never included in the BAM alignment file if they fail these guards
		bool     require_complete_match;  // Default = false   COMMAND-LINE OPTION
		uint32_t require_match_length;    // Default = 0 (OFF) COMMAND-LINE OPTION
		double   require_match_fraction;  // Default = 0.9     COMMAND-LINE OPTION
    //! ignore reads with this many or more mismatches (I+D+MM)
    int32_t  maximum_read_mismatches;     // Default = -1 (OFF)

    bool smalt;                       // Unused
		uint32_t max_smalt_diff;          // Unused
    
    //! Settings: Candidate Junction Prediction
    uint32_t preprocess_junction_min_indel_split_length;  // Default = 3
		uint32_t required_both_unique_length_per_side;        // Default = 5
		uint32_t required_one_unique_length_per_side;         // Default = ssaha2_seed_length = 13
    
    // Which candidate junctions do we test?
		uint32_t minimum_candidate_junction_pos_hash_score;   // Default = 2
		uint32_t maximum_junction_sequence_insertion_length;  // Default = 20
    uint32_t maximum_junction_sequence_overlap_length;    // Default = 20
		uint32_t minimum_candidate_junctions;                 // Default = 10
		uint32_t maximum_candidate_junctions;                 // Default = 5000
    double maximum_candidate_junction_length_factor;      // Default = 0.1
    bool penalize_negative_junction_overlap;              // Manually set. True for experimental treatment.
    
    //! Settings: Alignment Resolution
		bool add_split_junction_sides;                        // Default = true (possibly remove this option)  
    double junction_pos_hash_neg_log10_p_value_cutoff;    // Default = 3
    
    //! Settings: Mutation Identification
    
    //! ignore bases below this cutoff for RA evidence (still counted for deletions?)
    uint32_t base_quality_cutoff;                         // Default 3 COMMAND-LINE OPTION
        
    //! treated as read numbers if integers >= 1.0 and percentages of average coverage if > 0.0 and < 1.0
    double deletion_coverage_propagation_cutoff;          // Default = calculated COMMAND-LINE OPTION
    double deletion_coverage_seed_cutoff;                 // Default = 0;         COMMAND-LINE OPTION
    
    //! These are mutually exclusive settings (polymorphism prediction overrides mixed_base_prediction)
    bool polymorphism_prediction;                         // Default = false COMMAND-LINE OPTION
    //! Predict not only consensus genotype calls, but test mixed states between them.
    bool mixed_base_prediction;                           // Default = true
    
    //! Verbose output of bases encountered at each position
    bool mutation_identification_per_position_file;

    double mutation_log10_e_value_cutoff;                   // Default = 10
    double polymorphism_log10_e_value_cutoff;               // Default = mutation_log10_e_value_cutoff = 10
		double polymorphism_bias_p_value_cutoff;                // Default = 0.05
		double polymorphism_frequency_cutoff;                   // Default = 0.1 for mixed base | 0.0 for polymorphism
		int32_t polymorphism_minimum_new_coverage_each_strand; // Default = 1
		uint32_t polymorphism_reject_homopolymer_length;        // Default = 0 (OFF)
		bool no_indel_polymorphisms;                            // Default = false
		
		//! Settings: Copy Number Variation
    uint32_t copy_number_variation_tile_size;
    bool ignore_redundant_coverage;
    uint32_t periodicity_method;
    uint32_t periodicity_start;
    uint32_t periodicity_end;
    uint32_t periodicity_step;
    
    //! Settings: Output
    uint32_t maximum_reads_to_align;                      // Default = 100
    uint32_t max_rejected_polymorphisms_to_show;          // Default = 20
		uint32_t max_rejected_junctions_to_show;              // Default = 10
		bool hide_circular_genome_junctions;                  // Default = true (remove as option?)
    //! special output for Blount paper - not implemented in C++!
		bool lenski_format;                                   // Default = false (not implemented!)
    //! don't create any HTML evidence files
		bool no_evidence;                                     // Default = false (rarely used)
    //! show blue boxes rather than precentages in predictions
		bool shade_frequencies;                               // Default = false
    
    //! Settings: Experimental
    
    //! @GRC added in for gathering/analyzing breseq values
    bool add_metadata_to_gd;                              // Default = false COMMAND-LINE OPTION


    ////////////////////
    //! File Paths
    ////////////////////
    
    //! Paths populated from location of executable
		string bin_path;                  // absolute path to where this binary resides
    string program_data_path;         // path to where R scripts and images reside
		map<string, string> installed;    // hash of paths to programs that we call
    
    //! Paths: Sequence conversion
    string sequence_conversion_path;
    string sequence_conversion_done_file_name;

		string converted_fastq_file_name;
		string unwanted_fasta_file_name;
    string reference_trim_file_name;
		string sequence_conversion_summary_file_name;
    
    //! Paths: Read alignment
		string reference_alignment_path;
    string reference_alignment_done_file_name;
    
		string reference_hash_file_name;
		string reference_sam_file_name;
    
    //! Paths: Junction Prediction
    string candidate_junction_path;
    
    string preprocess_junction_done_file_name;
    string preprocess_junction_best_sam_file_name;
		string preprocess_junction_split_sam_file_name;
    
    string candidate_junction_done_file_name;
		string coverage_junction_best_bam_unsorted_file_name;
		string coverage_junction_best_bam_file_name;
		string coverage_junction_best_bam_prefix;
		string coverage_junction_distribution_file_name;
		string coverage_junction_plot_file_name;
		string coverage_junction_summary_file_name;
    string coverage_junction_error_count_summary_file_name;

    string coverage_junction_done_file_name;
    string candidate_junction_fasta_file_name;
    string candidate_junction_faidx_file_name;
		string candidate_junction_summary_file_name;
    
    //! Paths: Junction Alignment
    string candidate_junction_alignment_path;
    string candidate_junction_alignment_done_file_name;
    
    string candidate_junction_hash_file_name;
		string candidate_junction_sam_file_name;
    
    //! Paths: Alignment Resolution
    string alignment_resolution_path;
		string alignment_correction_done_file_name;
    
		string resolved_reference_sam_file_name;
		string resolved_junction_sam_file_name;
		string alignment_resolution_summary_file_name;
    string jc_genome_diff_file_name;
    
    //! Paths: BAM conversion
		string bam_path;
    string bam_done_file_name;

		string reference_bam_unsorted_file_name;
		string junction_bam_unsorted_file_name;
		string junction_bam_prefix;
		string junction_bam_file_name;

		//! Paths: Error Calibration
		string error_calibration_path;
		string error_counts_done_file_name;
		string error_rates_done_file_name;
    
		string error_counts_file_name;
		string error_rates_file_name;
		string coverage_file_name;
		string unique_only_coverage_distribution_file_name;
		string error_rates_summary_file_name;
		string error_rates_base_qual_error_prob_file_name;
		string plot_error_rates_r_script_file_name;
		string plot_error_rates_fit_r_script_file_name;
		string plot_error_rates_r_script_log_file_name;

		//! Paths: Mutation Identification
		string mutation_identification_path;
    
		string mutation_identification_done_file_name;    
		string complete_mutations_text_file_name;
		string complete_coverage_text_file_name;
		string genome_error_counts_file_name;
    string ra_mc_genome_diff_file_name;
    
    string polymorphism_statistics_done_file_name;
    string polymorphism_statistics_input_file_name;
		string polymorphism_statistics_output_file_name;
		string polymorphism_statistics_r_script_file_name;
		string polymorphism_statistics_r_script_log_file_name;
		string polymorphism_statistics_ra_mc_genome_diff_file_name;
    
		//! Paths: Copy Number Variation
    string copy_number_variation_path;
    string copy_number_variation_done_file_name;
    
    string tiled_complete_coverage_text_file_name;
    string ranges_text_file_name;
    string cnv_history_text_file_name;
    string smoothed_ranges_text_file_name;
    string final_cnv_text_file_name;
    string copy_number_variation_cn_genome_diff_file_name;
    
    string periodicity_table_file_name;
    string periodicity_done_file_name;
    
    //! Paths: Output
    string output_path;
		string output_done_file_name;

		string log_file_name;
		string index_html_file_name;
		string summary_html_file_name;
		string marginal_html_file_name;

    string local_evidence_path;
		string evidence_path;
    string evidence_genome_diff_file_name;
    string final_genome_diff_file_name;

		string local_coverage_plot_path;
		string coverage_plot_path;
    string coverage_plot_r_script_file_name;
		string overview_coverage_plot_file_name;
    
		string output_calibration_path;
		string unique_only_coverage_plot_file_name;
		string error_rates_plot_file_name;
    
		string breseq_small_graphic_from_file_name;
		string breseq_small_graphic_to_file_name;
    
    //! Paths: Data
    string data_path;
    
		string reference_bam_prefix;
		string reference_bam_file_name;
    string reference_fasta_file_name;
		string reference_faidx_file_name;
		string reference_gff3_file_name;
    string unmatched_read_file_name;
    
    //! Paths: Experimental
		string long_pairs_file_name;

    
    ////////////////////
    //! Methods
    ////////////////////
    
    // Set up defaults here
		Settings(const string& _base_output_path = "");
    
    // Constructor for default action
    Settings(int argc, char* argv[]);
    
    static void command_line_run_header();
    
		//! Utility function to substitute specific details into a generic file name
		static string file_name(const string& file_name_key, const string& substitute, const string& with)
		{
			string s(file_name_key);

			if (substitute.size() > 0)
			{
				size_t pos = s.find(substitute);
				if (pos != string::npos)
				{
					s.replace(pos, 1, with);
				}
			}

			return s;
		}
    
		//! Utility function to get relative path from two file locations for making HTML files
    static string relative_path(const string& full_path, const string& base_path) 
    {
      string s(full_path);
      if (s.substr(0, base_path.size()) == base_path)
      {
        s.erase(0,base_path.size());
      }
      
      if (s.at(0) == '/')
        s.erase(0,1);
      
      return s;
    }

		string ctool(string tool_name, bool allow_fail = false)
		{
      if (this->installed[tool_name].size() == 0)
      {
        cerr << "Executable '" << tool_name << "' not found in breseq bin path '" << this->bin_path << "'." << endl;
        if (!allow_fail) exit(-1);
        return "";
      }
			return tool_name;
		}
    
    static string time2string(const time_t& _time)
    {
      const struct tm * time_info = localtime(&_time);
      char s[1024];
      strftime(s, 1024, "%H:%M:%S %d %b %Y", time_info);
      return s;
    }
    
    static string elapsedtime2string(double _diff_time)
    {
      stringstream ss;
      
      double t = _diff_time; // in seconds
      uint32_t tm_yday = static_cast<uint32_t>(floor( t / (60*60*24)));
      t -= tm_yday * 60*60*24;
      uint32_t tm_hour = static_cast<uint32_t>(floor( t / (60*60)));
      t -= tm_hour * 60*60;
      uint32_t tm_min = static_cast<uint32_t>(floor( t / (60)));
      t -= tm_min * 60;
      uint32_t tm_sec = static_cast<uint32_t>(floor( t / (1)));      
      if (tm_yday > 0)
      {
        if (ss.str().length() > 0) ss << " ";
        ss << tm_yday;
        ss << ((tm_yday!=1)?" days":" day");
      }

      if (tm_hour > 0)
      {
        if (ss.str().length() > 0) ss << " ";
        ss << tm_hour;
        ss << ((tm_hour!=1)?" hours":" hour");
      }
      
      if (tm_min > 0)
      {
        if (ss.str().length() > 0) ss << " ";
        ss << tm_min;
        ss << ((tm_min!=1)?" minutes":" minute");
      }
      
      //if (tm_sec > 0)
      {
        if (ss.str().length() > 0) ss << " ";
        ss << tm_sec;
        ss << ((tm_sec!=1)?" seconds":" second");
      }

      return ss.str();
    }

    
    void record_start_time(const string& message)
    {
      ExecutionTime ex_time;
      time_t this_time = time(NULL);
      ex_time._time_start = this_time;
      ex_time._formatted_time_start = time2string(this_time);
      ex_time._time_end = 0;
      ex_time._formatted_time_end = "";
      ex_time._time_elapsed = 0;
      ex_time._formatted_time_elapsed = "";
      ex_time._message = message;

      this->execution_times.push_back(ex_time);
    }

		void record_end_time(const string& message)
		{
			uint32_t i = 0;
			while (i < this->execution_times.size())
			{
				if (this->execution_times[i]._message == message)
					break;
				i++;
			}

			if (i >= this->execution_times.size())
			{
				cout << "Did not find matching start time for:" << endl << message << endl;
				ExecutionTime blank;
				this->execution_times.push_back(blank);
			}

			ExecutionTime& ex_time = this->execution_times[i];
			ex_time._message = message;

			time_t this_time = time(NULL);
			ex_time._time_end = this_time;
			ex_time._formatted_time_end = time2string(this_time);

			//if we had a previous time, calculate elapsed
			if (i < this->execution_times.size())
			{
				ex_time._time_elapsed = difftime(ex_time._time_end, ex_time._time_start);
				ex_time._formatted_time_elapsed = elapsedtime2string(ex_time._time_elapsed);
			}
		}
    
 
		string base_name_to_read_file_name(string base_name)
		{
      return this->read_files.base_name_to_read_file_name(base_name);
		}

		bool do_step(string done_key, string message);
		void done_step(string done_key);
		void check_installed();
    void log(const string& message);
  private:

		void pre_option_initialize(int argc = 0, char* argv[] = NULL);
		void post_option_initialize();
    void init_installed();
  };

	class AnalyzeFastq : public Storable {
  public:
		uint32_t max_read_length;
		uint32_t num_reads;
		uint32_t min_quality_score;
		uint32_t max_quality_score;
		uint32_t num_bases;
		string original_qual_format;
    string quality_format;
		string converted_fastq_name;
    
    AnalyzeFastq() {};
    
    AnalyzeFastq(
                 uint32_t _max_read_length, 
                 uint32_t _num_reads, 
                 uint32_t _min_quality_score, 
                 uint32_t _max_quality_score, 
                 uint32_t _num_bases, 
                 const string& _original_qual_format, 
                 const string& _quality_format,
                 const string& _converted_fastq_name
                )
    : max_read_length(_max_read_length)
    , num_reads(_num_reads)
    , min_quality_score(_min_quality_score)
    , max_quality_score(_max_quality_score)
    , num_bases(_num_bases)
    , original_qual_format(_original_qual_format)
    , quality_format(_quality_format)
    , converted_fastq_name(_converted_fastq_name)
    { }
    
    void serialize(ofstream& f)
    {
      write_to_file(f, max_read_length);
      write_to_file(f, num_reads);
      write_to_file(f, min_quality_score);
      write_to_file(f, max_quality_score);
      write_to_file(f, num_bases);
      write_to_file(f, original_qual_format);
      write_to_file(f, quality_format);
      write_to_file(f, converted_fastq_name);
    }
    void deserialize(ifstream& f)
    {
      read_from_file(f, max_read_length);
      read_from_file(f, num_reads);
      read_from_file(f, min_quality_score);
      read_from_file(f, max_quality_score);
      read_from_file(f, num_bases);
      read_from_file(f, original_qual_format);
      read_from_file(f, quality_format);
      read_from_file(f, converted_fastq_name);
    }
	};

	class Summary : public Storable
	{
	public:

		class AlignmentResolution : public Storable
		{
    public:
      
      // @JEB TODO: these statistics are not implemented
      class ReadFile {
      public:
        uint32_t num_unmatched_reads;
        uint32_t num_unique_reads;
        uint32_t num_repeat_reads;
        
        ReadFile()
        : num_unmatched_reads(0)
        , num_unique_reads(0)
        , num_repeat_reads(0)
        {}
        
        void serialize(ofstream& f)
        {
          write_to_file(f, num_unmatched_reads);
          write_to_file(f, num_unique_reads);
          write_to_file(f, num_repeat_reads);
        }
        void deserialize(ifstream& f)
        {
          read_from_file(f, num_unmatched_reads);
          read_from_file(f, num_unique_reads);
          read_from_file(f, num_repeat_reads);
        }
        
      };
      storable_map<string,ReadFile> read_file;
      uint32_t total_unmatched_reads;
      uint32_t total_unique_reads;
      uint32_t total_repeat_reads;
      
      // these statistics are implemented
      
      //! map by reference seq_id of number of bases from possible overlap
      //! that will be accepted
      map<string,int32_t> distance_cutoffs; 
      //! map by reference seq_id, then list by non-overlap distance possible
      //! of minimum pos hash score needed for accepting a junction
      storable_map<string, storable_vector<int32_t> > pos_hash_cutoffs;   
      map<int32_t, int32_t> observed_pos_hash_score_distribution;
      map<int32_t, int32_t> accepted_pos_hash_score_distribution;

			void serialize(ofstream& f)
			{
        read_file.serialize(f);
        write_to_file(f, distance_cutoffs);
				pos_hash_cutoffs.serialize(f);
				write_to_file(f, observed_pos_hash_score_distribution);
				write_to_file(f, accepted_pos_hash_score_distribution);
			}
			void deserialize(ifstream& f)
			{
        read_file.deserialize(f);
        read_from_file(f, distance_cutoffs);
        pos_hash_cutoffs.deserialize(f);
				read_from_file(f, observed_pos_hash_score_distribution);
				read_from_file(f, accepted_pos_hash_score_distribution);
			}

		} alignment_resolution;

		storable_map<string, Coverage> preprocess_coverage;
		storable_map<string, Coverage> unique_coverage;

		class CandidateJunctionSummaryData : public Storable
    {
    public:
			struct Total
			{
				int32_t number;
				int32_t length;
			} total;

			struct Accepted
			{
				int32_t number;
				int32_t length;
				int32_t pos_hash_score_cutoff;
			} accepted;

			map<int32_t, int32_t> pos_hash_score_distribution;

      void serialize(ofstream& f)
      {
        write_to_file(f, total);
        write_to_file(f, accepted);
				write_to_file(f, pos_hash_score_distribution);
      }
      
      void deserialize(ifstream& f)
      {
        read_from_file(f, total);
        read_from_file(f, accepted);
				read_from_file(f, pos_hash_score_distribution);
      }

		} candidate_junction;

		class SequenceConversion : public Storable
		{
    public:
			float avg_read_length;
			uint32_t max_qual;
			uint32_t num_reads;
			uint32_t num_bases;
			map<string, string> converted_fastq_name;
			storable_map<string, AnalyzeFastq> reads;
			uint32_t total_reference_sequence_length;
			uint32_t max_read_length;

			void serialize(ofstream& f)
			{
				write_to_file(f, avg_read_length);
        write_to_file(f, max_qual);
				write_to_file(f, num_reads);
				write_to_file(f, num_bases);
				write_to_file(f, converted_fastq_name);
        reads.serialize(f);
        write_to_file(f, total_reference_sequence_length);
				write_to_file(f, max_read_length);
			}
			void deserialize(ifstream& f)
			{
        read_from_file(f, avg_read_length);
        read_from_file(f, max_qual);
				read_from_file(f, num_reads);
				read_from_file(f, num_bases);
				read_from_file(f, converted_fastq_name);
        reads.deserialize(f);
        read_from_file(f, total_reference_sequence_length);
				read_from_file(f, max_read_length);
			}

		} sequence_conversion;
    
    class ErrorCount : public Storable
    {
    public:
      double no_pos_hash_per_position_pr;
      
      void serialize(ofstream& f)
      {
        write_to_file(f, no_pos_hash_per_position_pr);
      }
      void deserialize(ifstream& f)
      {
        read_from_file(f, no_pos_hash_per_position_pr);
      }
    };
    
    storable_map<string, ErrorCount> preprocess_error_count;


    // Overall functions for all of summary

		void serialize(ofstream& f)
		{
      sequence_conversion.serialize(f);
      candidate_junction.serialize(f);
      alignment_resolution.serialize(f);
      preprocess_coverage.serialize(f);
      unique_coverage.serialize(f);
      preprocess_error_count.serialize(f);
    }
    
		void deserialize(ifstream& f)
		{
      sequence_conversion.deserialize(f);
      candidate_junction.deserialize(f);
      alignment_resolution.deserialize(f);
      preprocess_coverage.deserialize(f);
      unique_coverage.deserialize(f);
      preprocess_error_count.deserialize(f);
		}
	};
  
} // breseq namespace

#endif
