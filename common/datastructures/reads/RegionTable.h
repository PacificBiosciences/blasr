#ifndef DATASTRUCTURES_READS_REGION_TABLE_H_
#define DATASTRUCTURES_READS_REGION_TABLE_H_

#include <assert.h>
#include <string>
#include <vector>
#include "../../Enumerations.h"
using namespace std;

/*
  It appears that this class is being passed to something in the STL and that it is being compared to an int.
  That comparison fails b/c there is no explicit cast-to-int for a RegionAnnotation.
  
  It looks like either you intend for some member function's return value to be passed to the STL algorithm,
  or you want to write RegionAnnotation::int()const {} to allow the implicit cast to work
 */
class RegionAnnotation {
 public:
	typedef enum T_AnnotationRow {HoleNumber, RegionType, RegionStart, RegionEnd, RegionScore} AnnotationRow;
	static const int NCOLS=5;
	int row[NCOLS];
	int operator<(const RegionAnnotation &rhs) const {
		return row[0] < rhs.row[0];
	}

	int operator<(int holeNumber) {
		return row[0] < holeNumber;
	}

	RegionAnnotation& operator=(const RegionAnnotation &rhs) {
		memcpy(row, rhs.row, sizeof(int)*NCOLS);
    return *this;
	}
  int GetHoleNumber() {
    return row[HoleNumber];
  }

  int SetHoleNumber(int holeNumber) {
    row[HoleNumber] = holeNumber;
  }

  int GetType() {
    return row[RegionType];
  }
  
  int SetType(int regionType) {
    row[RegionType] = regionType;
  }

  int GetStart() {
    return row[RegionStart];
  }

  void SetStart(int start) {
    row[RegionStart] = start;
  }
  int GetEnd() {
    return row[RegionEnd];
  }
  
  void SetEnd(int end) {
    row[RegionEnd] = end;
  }
  
  int GetScore() {
    return row[RegionScore];
  }
  
  void SetScore(int score) {
    row[RegionScore] = score;
  }
  
};

class RegionTable {
 public:
	vector<RegionAnnotation> table;
	vector<string> columnNames;
	vector<string> regionTypes;
	vector<string> regionDescriptions;
	vector<string> regionSources;
	vector<RegionType> regionTypeEnums;

	int LookupRegionsByHoleNumber(int holeNumber, int &low, int &high) {
		vector<RegionAnnotation>::iterator lowIt, highIt;
		lowIt  = std::lower_bound(table.begin(), table.end(), holeNumber);
		highIt = std::lower_bound(table.begin(), table.end(), holeNumber+1);
		low =  lowIt - table.begin();
		high = highIt - table.begin();
		return high-low;
	}

  //
  // Define a bunch of accessor functions.
  //

  //
  // Different region tables have different ways of encoding regions.
  // This maps from the way they are encoded in the rgn table to a
  // standard encoding.
  //

  RegionType GetType(int regionIndex) {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    return (RegionType) regionTypeEnums[table[regionIndex].GetType()];
  }
  
  int GetStart(int regionIndex) {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    return table[regionIndex].GetStart();
  }
  
  void SetStart(int regionIndex, int start) {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    table[regionIndex].SetStart(start);
  }
  
  int GetEnd(int regionIndex) {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    return table[regionIndex].GetEnd();
  }

  void SetEnd(int regionIndex, int end) {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    table[regionIndex].SetEnd(end);
  }

  int GetHoleNumber(int regionIndex) {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    return table[regionIndex].GetHoleNumber();
  }

  int SetHoleNumber(int regionIndex, int holeNumber) {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    return table[regionIndex].SetHoleNumber(holeNumber);
  }

  int GetScore(int regionIndex) {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    return table[regionIndex].row[RegionAnnotation::RegionScore];
  }

  int SetScore(int regionIndex, int score) {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    return table[regionIndex].row[RegionAnnotation::RegionScore] = score;
  }

	void SortTableByHoleNumber() {
		std::sort(table.begin(), table.end());
	}

	void Reset() {
		table.clear();
		columnNames.clear();
		regionTypes.clear();
		regionDescriptions.clear();
		regionSources.clear();
		regionTypeEnums.clear();
	}

  void CreateDefaultAttributes() {
    columnNames.clear();
    columnNames.push_back("HoleNumber");
    columnNames.push_back("Region type index");
    columnNames.push_back("Region start in bases");
    columnNames.push_back("Region end in bases");
    columnNames.push_back("Region score");

    regionTypes.push_back("Adapter");
    regionTypes.push_back("Insert");
    regionTypes.push_back("HQRegion");
    
    regionDescriptions.push_back("Adapter Hit");
    regionDescriptions.push_back("Insert Region");
    regionDescriptions.push_back("High Quality bases region. Score is 1000 * "
                                 "predicted accuracy, where predicted accuary is 0 to 1.0"); 
    
    regionSources.push_back("AdapterFinding");
    regionSources.push_back("AdapterFinding");
    regionSources.push_back("PulseToBase Region classifer");

    regionTypeEnums.push_back(Adapter);
    regionTypeEnums.push_back(Insert);
    regionTypeEnums.push_back(HQRegion);
  }
};


#endif
