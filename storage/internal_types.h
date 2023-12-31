
#pragma once
#pragma GCC diagnostic ignored "-Wcast-qual"

//#define UNUSED(x) (void)(x)
#include <unistd.h>
#include <bitset>
#include <climits>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "my_inttypes.h"

#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_vector.h"

//===------------nvm--------------------------------------------------------===//
#include <unistd.h>
#include <iostream>
#include <libpmemobj++/make_persistent.hpp>
#include <libpmemobj++/p.hpp>
#include <libpmemobj++/persistent_ptr.hpp>
#include <libpmemobj++/transaction.hpp>
#include <libpmemobj.h>
#include <libpmemobj++/make_persistent_atomic.hpp>
#include <libpmemobj++/pool.hpp>
#include <libpmemobj++/container/concurrent_hash_map.hpp>
#include <libpmemobj++/container/basic_string.hpp>
#include <libpmemobj++/container/vector.hpp>


#define COLUMN_STORE_LAYOUT_OID 1
#define INVALID_TYPE_ID 0
#define ROW_STORE_LAYOUT_OID 0
//typedef unsigned int uint32_t;
//typedef signed __int64       int64_t;
//typedef unsigned __int64     uint64_t;

/* column_map_type used to store the layout of a tile_group
 * <column offset> to <tile offset, tile column offset>
 */
typedef std::map<uint, std::pair<uint, uint>> column_map_type;
static const uint INVALID_OID = std::numeric_limits<uint>::max();

//==================nvm===============
#define TILENVMSIZE ((size_t) (1024 * 1024 * 64)) /* 64 MB */
struct Tile_nvm {
  pmem::obj::persistent_ptr<Tile_nvm> next;
  pmem::obj::p<size_t> offset;//当前插入了多少条记录
//  char *data;
//  pmem::obj::persistent_ptr<uchar*> *data;
  pmem::obj::persistent_ptr<char[]> data;
//    char data[500000];

};
using block_tile_type = pmem::obj::concurrent_hash_map<pmem::obj::p<uint64_t>,
                                                       pmem::obj::persistent_ptr<Tile_nvm>>;


//extern Tile_nvm tile_nvm;
//============nvm===============
struct root {
  pmem::obj::p<size_t> num;
  pmem::obj::persistent_ptr<Tile_nvm> tile_nvms;
  pmem::obj::persistent_ptr<block_tile_type> block_tile;
};

//1,000,000* 200
static const uint64_t INDEX_LOG_TILE_LEN = 1000000*200;
static const uint64_t INDEX_LOG_TILE_COUNT = 1000;
static const uint64_t INDEX_LOG_COUNT = 1000;

using chpt_map_type = pmem::obj::concurrent_hash_map<pmem::obj::p<uint64_t>,
                                                pmem::obj::p<uint64_t>>;


struct chpt_tuple {
  pmem::obj::p<uint64_t> key;
  pmem::obj::p<uint64_t> point;
};
using vector_type = pmem::obj::vector<chpt_tuple>;

struct chpt_table {
  pmem::obj::p<uint32_t> num;
  pmem::obj::p<uint32_t> next_chptcopy = 1;
  pmem::obj::persistent_ptr<chpt_tuple[]>  chpt_table1;
  pmem::obj::persistent_ptr<chpt_map_type> chpt_map1;
  pmem::obj::persistent_ptr<chpt_map_type> chpt_map2;
};

struct Index_log {

  pmem::obj::p<uint64_t> key;//四个16位的数字组成.
  pmem::obj::p<uint32_t> block;//itempoint -1则为delete
  pmem::obj::p<uint32_t> offset;//itempoint -1则为delete
  pmem::obj::p<uint64_t> txid;
  pmem::obj::p<uint8_t> tid;//索引id 或者tableid
};

struct Index_log_tile {
  pmem::obj::p<uint32_t> num;
  pmem::obj::persistent_ptr<Index_log_tile> next;
//  pmem::obj::persistent_ptr<index_log[]> data;
  pmem::obj::persistent_ptr<Index_log[]> data;// 1000 Index_log/data

};

//using hashmap_type = pmem::obj::concurrent_hash_map<pmem::obj::p<uint64_t>, pmem::obj::p<uint64_t>>;

using hashmap_type = pmem::obj::concurrent_hash_map<pmem::obj::p<uint64_t>,
                                                    pmem::obj::persistent_ptr<Index_log_tile>>;


struct root_index_log {
  pmem::obj::p<size_t> num;
  pmem::obj::persistent_ptr<Index_log_tile> curr_index_log_tile;
  pmem::obj::persistent_ptr<hashmap_type> tx_idx_log_tile;
};

struct Txid_idx_log_tile {
//  pmem::obj::persistent_ptr<std::unordered_map<uint64_t ,uint64_t>> num;
  pmem::obj::persistent_ptr<hashmap_type> root_tile;
};

//extern root root_nvm;
// Root structure
struct my_root {
  size_t len;
  char buf[20];
};

//using col_type_nvm = pmem::obj::concurrent_hash_map<pmem::obj::persistent_ptr<char[]>,
//                                                    pmem::obj::p<uint32_t>>;

struct index_nvm {
  pmem::obj::persistent_ptr<char[]> index_name;
  pmem::obj::p<uint32_t> index_type;//1 pri 2 unique
  pmem::obj::p<uint32_t> index_len;
  pmem::obj::p<uint32_t> col_no;
  pmem::obj::persistent_ptr<uint32_t []> index_cols;
  pmem::obj::persistent_ptr<index_nvm> next;
  pmem::obj::p<uint32_t> index_oid;


};

struct col_nvm {
  pmem::obj::persistent_ptr<char[]> col_name;
  pmem::obj::p<uint32_t> col_type;
};

struct table_nvm {
  pmem::obj::persistent_ptr<char[]> name;
  pmem::obj::persistent_ptr<table_nvm> next;
  pmem::obj::persistent_ptr<col_nvm[]> col_type;//col<name, type>
  //index name, cols(1,2,3 )
  pmem::obj::p<uint32_t> col_no = 0;
  pmem::obj::p<uint16_t> move_pos = 0;
  pmem::obj::p<uint16_t> tableId = 0;
  pmem::obj::persistent_ptr<index_nvm> root_index_nvm;
//  pmem::obj::persistent_ptr<block_tile_type> block_tile;
//  pmem::obj::persistent_ptr<Tile_nvm> tile_nvms;
  //1: uint 2: int64 3: double 4: varchar
};



struct database_nvm {
  pmem::obj::persistent_ptr<char[]> name;
  pmem::obj::p<uint32_t> table_num = 0;
  pmem::obj::persistent_ptr<table_nvm> root_table;
  //1: uint 2: int64 3: double 4: varchar
};



// For transaction id

typedef uint64_t txn_id_t;


// For commit id

typedef uint64_t cid_t;



// For epoch id

typedef uint64_t eid_t;
//typedef uint64_t txn_id_t;
static const uint START_OID = 0;
static const txn_id_t INVALID_TXN_ID = 0;
static const cid_t INVALID_CID = 0;



//static const txn_id_t INITIAL_TXN_ID = 1;
static const txn_id_t INITIAL_TXN_ID = 1;

//int DEFAULT_TUPLES_PER_TILEGROUP = 1000;
//extern int DEFAULT_TUPLES_PER_TILEGROUP;

static const cid_t MAX_CID = std::numeric_limits<cid_t>::max();

//static const txn_id_t MAX_TXN_ID = std::numeric_limits<txn_id_t>::max();
static const txn_id_t MAX_TXN_ID = std::numeric_limits<txn_id_t>::max();

enum class LayoutType {
    INVALID = INVALID_TYPE_ID,
    ROW = 1,    /* Pure row layout */
    COLUMN = 2, /* Pure column layout */
    HYBRID = 3  /* Hybrid layout */
};

enum class BackendType {
    INVALID = INVALID_TYPE_ID,  // invalid backend type
    MM = 1,                     // on volatile memory
    NVM = 2,                    // on non-volatile memory
    SSD = 3,                    // on ssd
    HDD = 4                     // on hdd
};

enum class ConstraintType {
    INVALID = INVALID_TYPE_ID,  // invalid
    CHECK = 1,                  // check
    PRIMARY = 2,                // primary key
    UNIQUE = 3,                 // unique
    FOREIGN = 4,                // foreign key
    EXCLUSION = 5               // foreign key
};

enum class FKConstrActionType {
    INVALID = INVALID_TYPE_ID,  // invalid
    NOACTION = 1,
    RESTRICT = 2,
    CASCADE = 3,
    SETNULL = 4,
    SETDEFAULT = 5
};

enum class StartArg {
  INVALID = INVALID_TYPE_ID,
  NO_TABLE = 1,
  HAS_TABLE = 2,
  HAS_INDEX = 3,
};


enum class ExpressionType {
    INVALID = INVALID_TYPE_ID,

    // -----------------------------
    // Arithmetic Operators
    // Implicit Numeric Casting: Trying to implement SQL-92.
    // Implicit Character Casting: Trying to implement SQL-92, but not easy...
    // Anyway, use explicit OPERATOR_CAST if you could.
    // -----------------------------

    // left + right (both must be number. implicitly casted)
    OPERATOR_PLUS = 1,
    // left - right (both must be number. implicitly casted)
    OPERATOR_MINUS = 2,
    // left * right (both must be number. implicitly casted)
    OPERATOR_MULTIPLY = 3,
    // left / right (both must be number. implicitly casted)
    OPERATOR_DIVIDE = 4,
    // left || right (both must be char/varchar)
    OPERATOR_CONCAT = 5,
    // left % right (both must be integer)
    OPERATOR_MOD = 6,
    // explicitly cast left as right (right is integer in ValueType enum)
    OPERATOR_CAST = 7,
    // logical not operator
    OPERATOR_NOT = 8,
    // is null operator
    OPERATOR_IS_NULL = 21,
    // is not null operator
    OPERATOR_IS_NOT_NULL = 22,
    // exists test.
    OPERATOR_EXISTS = 18,
    OPERATOR_UNARY_MINUS = 60,

    // -----------------------------
    // Comparison Operators
    // -----------------------------
    // equal operator between left and right
    COMPARE_EQUAL = 10,
    // inequal operator between left and right
    COMPARE_NOTEQUAL = 11,
    // less than operator between left and right
    COMPARE_LESSTHAN = 12,
    // greater than operator between left and right
    COMPARE_GREATERTHAN = 13,
    // less than equal operator between left and right
    COMPARE_LESSTHANOREQUALTO = 14,
    // greater than equal operator between left and right
    COMPARE_GREATERTHANOREQUALTO = 15,
    // LIKE operator (left LIKE right). Both children must be string.
    COMPARE_LIKE = 16,
    // NOT LIKE operator (left NOT LIKE right). Both children must be string.
    COMPARE_NOTLIKE = 17,
    // IN operator [left IN (right1, right2, ...)]
    COMPARE_IN = 19,
    // IS DISTINCT FROM operator
    COMPARE_DISTINCT_FROM = 20,

    // -----------------------------
    // Conjunction Operators
    // -----------------------------
    CONJUNCTION_AND = 30,
    CONJUNCTION_OR = 31,

    // -----------------------------
    // Values
    // -----------------------------
    VALUE_CONSTANT = 40,
    VALUE_PARAMETER = 41,
    VALUE_TUPLE = 42,
    VALUE_TUPLE_ADDRESS = 43,
    VALUE_NULL = 44,
    VALUE_VECTOR = 45,
    VALUE_SCALAR = 46,

    // -----------------------------
    // Aggregates
    // -----------------------------
    AGGREGATE_COUNT = 50,
    AGGREGATE_COUNT_STAR = 51,
    AGGREGATE_SUM = 52,
    AGGREGATE_MIN = 53,
    AGGREGATE_MAX = 54,
    AGGREGATE_AVG = 55,

    // -----------------------------
    // Functions
    // -----------------------------
    FUNCTION = 100,

    // -----------------------------
    // Internals added for Elastic
    // -----------------------------
    HASH_RANGE = 200,

    // -----------------------------
    // Operators
    // -----------------------------
    OPERATOR_CASE_EXPR = 302,
    OPERATOR_NULLIF = 304,
    OPERATOR_COALESCE = 305,

    // -----------------------------
    // Subquery IN/EXISTS
    // -----------------------------
    ROW_SUBQUERY = 400,
    SELECT_SUBQUERY = 401,

    // -----------------------------
    // Parser
    // -----------------------------
    STAR = 500,
    PLACEHOLDER = 501,
    COLUMN_REF = 502,
    FUNCTION_REF = 503,
    TABLE_REF = 504,

    // -----------------------------
    // Miscellaneous
    // -----------------------------
    CAST = 600
};

// Note that we have some duplicate DatePartTypes with the 's' suffix
// They have to have the same value in order for it to work.
enum class DatePartType : uint32_t {
    INVALID = INVALID_TYPE_ID,
    CENTURY = 1,
    DAY = 2,
    DAYS = 2,
    DECADE = 3,
    DECADES = 3,
    DOW = 4,
    DOY = 5,
    HOUR = 7,
    HOURS = 7,
    MICROSECOND = 10,
    MICROSECONDS = 10,
    MILLENNIUM = 11,
    MILLISECOND = 12,
    MILLISECONDS = 12,
    MINUTE = 13,
    MINUTES = 13,
    MONTH = 14,
    MONTHS = 14,
    QUARTER = 15,
    QUARTERS = 15,
    SECOND = 16,
    SECONDS = 16,
    WEEK = 20,
    WEEKS = 20,
    YEAR = 21,
    YEARS = 21,
};

std::string FKConstrActionTypeToString(FKConstrActionType type);
std::string ConstraintTypeToString(ConstraintType type);
std::string TypeIdToString(TypeId type);
TypeId StringToTypeId(const std::string &str);

//#endif



class ItemPointer;
struct ItemPointerHasher;
class ItemPointerComparator;

// For all of the enums defined in this header, we will
// use this value to indicate that it is an invalid value
// I don't think it matters whether this is 0 or -1
#define INVALID_TYPE_ID 0

//===--------------------------------------------------------------------===//
// NULL-related Constants
//===--------------------------------------------------------------------===//

#define VALUE_COMPARE_LESSTHAN -1
#define VALUE_COMPARE_EQUAL 0
#define VALUE_COMPARE_GREATERTHAN 1
#define VALUE_COMPARE_INVALID -2
#define VALUE_COMPARE_NO_EQUAL \
  -3  // assigned when comparing array list and no element matching.

#define INVALID_RATIO -1

#define DEFAULT_DB_NAME "default_database"

extern int DEFAULT_TUPLES_PER_TILEGROUP;
extern int TEST_TUPLES_PER_TILEGROUP;

//===--------------------------------------------------------------------===//
// Other Constants
//===--------------------------------------------------------------------===//

#define VARCHAR_LENGTH_SHORT 16
#define VARCHAR_LENGTH_MID 256
#define VARCHAR_LENGTH_LONG 4096

//===--------------------------------------------------------------------===//
// Comparison Result
// This is used when comparing values with each other
//===--------------------------------------------------------------------===//
enum class CmpBool {
  CmpFalse = 0,
  CmpTrue = 1,
  NULL_ = 2  // Note the underscore suffix
};

//===--------------------------------------------------------------------===//
// Postgres Value Types
// This file defines all the types that we will support
// We do not allow for user-defined types, nor do we try to do anything dynamic.
//
// For more information, see 'pg_type.h' in Postgres
// https://github.com/postgres/postgres/blob/master/src/include/catalog/pg_type.h#L273
//===--------------------------------------------------------------------===//

enum class PostgresValueType {
  INVALID = INVALID_TYPE_ID,
  BOOLEAN = 16,
  TINYINT = 16,  // BOOLEAN is an alias for TINYINT
  SMALLINT = 21,
  INTEGER = 23,
  VARBINARY = 17,
  BIGINT = 20,
  REAL = 700,
  DOUBLE = 701,
  TEXT = 25,
  BPCHAR = 1042,
  BPCHAR2 = 1014,
  VARCHAR = 1015,
  VARCHAR2 = 1043,
  DATE = 1082,
  TIMESTAMPS = 1114,
  TIMESTAMPS2 = 1184,
  TEXT_ARRAY = 1009,     // TEXTARRAYOID in postgres code
  INT2_ARRAY = 1005,     // INT2ARRAYOID in postgres code
  INT4_ARRAY = 1007,     // INT4ARRAYOID in postgres code
  OID_ARRAY = 1028,      // OIDARRAYOID in postgres code
  FLOADT4_ARRAY = 1021,  // FLOADT4ARRAYOID in postgres code
  DECIMAL = 1700
};
std::string PostgresValueTypeToString(PostgresValueType type);
PostgresValueType StringToPostgresValueType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const PostgresValueType &type);



// When short_str is true, return a short version of ExpressionType string
// For example, + instead of Operator_Plus. It's used to generate the
// expression name
std::string ExpressionTypeToString(ExpressionType type, bool short_str = false);
ExpressionType StringToExpressionType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const ExpressionType &type);
ExpressionType ParserExpressionNameToExpressionType(const std::string &str);

std::string DatePartTypeToString(DatePartType type);
DatePartType StringToDatePartType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const DatePartType &type);

// PAVLO: 2017-01-18
// I removed these DatePartTypes because I don't think that
// it's something we can easily support right now. Things get
// weird with timezones
//   EPOCH = 6,
//   ISODOW = 8,
//   ISOYEAR = 9,
//   TIMEZONE = 17,
//   TIMEZONE_HOUR = 18,
//   TIMEZONE_HOURS = 18,
//   TIMEZONE_MINUTE = 19,
//   TIMEZONE_MINUTES = 19,

//===--------------------------------------------------------------------===//
// Network Message Types
//===--------------------------------------------------------------------===//

enum class NetworkMessageType : unsigned char {
  // Important: The character '0' is treated as a null message
  // That means we cannot have an invalid type
  NULL_COMMAND = '0',

  // Responses
  PARSE_COMPLETE = '1',
  BIND_COMPLETE = '2',
  CLOSE_COMPLETE = '3',
  COMMAND_COMPLETE = 'C',
  PARAMETER_STATUS = 'S',
  AUTHENTICATION_REQUEST = 'R',
  ERROR_RESPONSE = 'E',
  EMPTY_QUERY_RESPONSE = 'I',
  NO_DATA_RESPONSE = 'n',
  READY_FOR_QUERY = 'Z',
  ROW_DESCRIPTION = 'T',
  DATA_ROW = 'D',
  // Errors
  HUMAN_READABLE_ERROR = 'M',
  SQLSTATE_CODE_ERROR = 'C',
  // Commands
  EXECUTE_COMMAND = 'E',
  SYNC_COMMAND = 'S',
  TERMINATE_COMMAND = 'X',
  DESCRIBE_COMMAND = 'D',
  BIND_COMMAND = 'B',
  PARSE_COMMAND = 'P',
  SIMPLE_QUERY_COMMAND = 'Q',
  CLOSE_COMMAND = 'C',
  // SSL willingness
  SSL_YES = 'S',
  SSL_NO = 'N',
};

enum class NetworkTransactionStateType : unsigned char {
  INVALID = static_cast<unsigned char>(INVALID_TYPE_ID),
  IDLE = 'I',
  BLOCK = 'T',
  FAIL = 'E',
};

enum class SqlStateErrorCode {
  SERIALIZATION_ERROR = '1',
};

//===--------------------------------------------------------------------===//
// Concurrency Control Types
//===--------------------------------------------------------------------===//

enum class ProtocolType {
  INVALID = INVALID_TYPE_ID,
  TIMESTAMP_ORDERING = 1  // timestamp ordering
};
std::string ProtocolTypeToString(ProtocolType type);
ProtocolType StringToProtocolType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const ProtocolType &type);

//===--------------------------------------------------------------------===//
// Epoch Types
//===--------------------------------------------------------------------===//

enum class EpochType {
  INVALID = INVALID_TYPE_ID,
  DECENTRALIZED_EPOCH = 1  // decentralized epoch manager
};
std::string EpochTypeToString(EpochType type);
EpochType StringToEpochType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const EpochType &type);

enum class TimestampType1 {
  INVALID = INVALID_TYPE_ID,
  SNAPSHOT_READ = 1,
  READ = 2,
  COMMIT = 3
};
std::string TimestampTypeToString(TimestampType1 type);
TimestampType1 StringToTimestampType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const TimestampType1 &type);

//===--------------------------------------------------------------------===//
// Visibility Types
//===--------------------------------------------------------------------===//

enum class VisibilityType {
  INVALID = INVALID_TYPE_ID,
  INVISIBLE = 1,
  DELETED = 2,
  OK = 3
};
std::string VisibilityTypeToString(VisibilityType type);
VisibilityType StringToVisibilityType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const VisibilityType &type);

enum class VisibilityIdType {
  INVALID = INVALID_TYPE_ID,
  READ_ID = 1,
  COMMIT_ID = 2
};
std::string VisibilityIdTypeToString(VisibilityIdType type);
VisibilityIdType StringToVisibilityIdType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const VisibilityIdType &type);

//===--------------------------------------------------------------------===//
// Isolation Levels
//===--------------------------------------------------------------------===//

enum class IsolationLevelType {
  INVALID = INVALID_TYPE_ID,
  SERIALIZABLE = 1,      // serializable
  SNAPSHOT = 2,          // snapshot isolation
  REPEATABLE_READS = 3,  // repeatable reads
  READ_COMMITTED = 4,    // read committed
};
std::string IsolationLevelTypeToString(IsolationLevelType type);
IsolationLevelType StringToIsolationLevelType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const IsolationLevelType &type);

//===--------------------------------------------------------------------===//
// Conflict Avoidance types
//===--------------------------------------------------------------------===//

enum class ConflictAvoidanceType {
  INVALID = INVALID_TYPE_ID,
  WAIT = 1,   // wait-based
  ABORT = 2,  // abort-based
};
std::string ConflictAvoidanceTypeToString(ConflictAvoidanceType type);
ConflictAvoidanceType StringToConflictAvoidanceType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const ConflictAvoidanceType &type);

//===--------------------------------------------------------------------===//
// Garbage Collection Types
//===--------------------------------------------------------------------===//

enum class GarbageCollectionType {
  INVALID = INVALID_TYPE_ID,
  OFF = 1,  // turn off GC
  ON = 2    // turn on GC
};
std::string GarbageCollectionTypeToString(GarbageCollectionType type);
GarbageCollectionType StringToGarbageCollectionType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const GarbageCollectionType &type);

//===--------------------------------------------------------------------===//
// Backend Types
//===--------------------------------------------------------------------===//


std::string BackendTypeToString(BackendType type);
BackendType StringToBackendType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const BackendType &type);

//===--------------------------------------------------------------------===//
// Index Types
//===--------------------------------------------------------------------===//

enum class IndexType {
  INVALID = INVALID_TYPE_ID,  // invalid index type
  BWTREE = 1,                 // bwtree
  HASH = 2,                   // hash
  SKIPLIST = 3,               // skiplist
  ART = 4,                    // ART
  //jy
  MASSTREE = 5,               //MASSTREE
};
std::string IndexTypeToString(IndexType type);
IndexType StringToIndexType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const IndexType &type);

enum class IndexConstraintType {
  // invalid index constraint type
  INVALID = INVALID_TYPE_ID,
  // default type - not used to enforce constraints
  DEFAULT = 1,
  // used to enforce primary key constraint
  PRIMARY_KEY = 2,
  // used for unique constraint
  UNIQUE = 3
};
std::string IndexConstraintTypeToString(IndexConstraintType type);
IndexConstraintType StringToIndexConstraintType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const IndexConstraintType &type);

//===--------------------------------------------------------------------===//
// Hybrid Scan Types
//===--------------------------------------------------------------------===//

enum class HybridScanType {
  INVALID = INVALID_TYPE_ID,
  SEQUENTIAL = 1,
  INDEX = 2,
  HYBRID = 3
};
std::string HybridScanTypeToString(HybridScanType type);
HybridScanType StringToHybridScanType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const HybridScanType &type);

//===--------------------------------------------------------------------===//
// Parse Node Types
//===--------------------------------------------------------------------===//

enum class ParseNodeType {
  INVALID = INVALID_TYPE_ID,  // invalid parse node type

  // Scan Nodes
  SCAN = 10,

  // DDL Nodes
  CREATE = 20,
  DROP = 21,

  // Mutator Nodes
  UPDATE = 30,
  INSERT = 31,
  DELETE = 32,

  // Prepared Nodes
  PREPARE = 40,
  EXECUTE = 41,

  // Select Nodes
  SELECT = 50,
  JOIN_EXPR = 51,  // a join tree
  TABLE = 52,      // a single table

  // Test
  MOCK = 80
};
std::string ParseNodeTypeToString(ParseNodeType type);
ParseNodeType StringToParseNodeType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const ParseNodeType &type);

//===--------------------------------------------------------------------===//
// Plan Node Types
//===--------------------------------------------------------------------===//

enum class PlanNodeType {
  INVALID = INVALID_TYPE_ID,  // invalid plan node type

  // Scan Nodes
  SEQSCAN = 10,
  INDEXSCAN = 11,
  CSVSCAN = 12,

  // Join Nodes
  NESTLOOP = 20,
  NESTLOOPINDEX = 21,
  MERGEJOIN = 22,
  HASHJOIN = 23,

  // Mutator Nodes
  UPDATE = 30,
  INSERT = 31,
  DELETE = 32,

  // DDL Nodes
  DROP = 33,
  CREATE = 34,
  POPULATE_INDEX = 35,
  ANALYZE = 36,

  // Communication Nodes
  SEND = 40,
  RECEIVE = 41,
  PRINT = 42,

  // Algebra Nodes
  AGGREGATE = 50,
  UNION = 52,
  ORDERBY = 53,
  PROJECTION = 54,
  MATERIALIZE = 55,
  LIMIT = 56,
  DISTINCT = 57,
  SETOP = 58,   // set operation
  APPEND = 59,  // append
  AGGREGATE_V2 = 61,
  HASH = 62,

  // Utility
  RESULT = 70,
  EXPORT_EXTERNAL_FILE = 71,
  CREATE_FUNC = 72,

  // Test
  MOCK = 80
};
std::string PlanNodeTypeToString(PlanNodeType type);
PlanNodeType StringToPlanNodeType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const PlanNodeType &type);

//===--------------------------------------------------------------------===//
// Create Types
//===--------------------------------------------------------------------===//

enum class CreateType {
  INVALID = INVALID_TYPE_ID,  // invalid create type
  DB = 1,                     // db create type
  TABLE = 2,                  // table create type
  INDEX = 3,                  // index create type
  CONSTRAINT = 4,             // constraint create type
  TRIGGER = 5,                // trigger create type
  SCHEMA = 6,                 // schema create type
};
std::string CreateTypeToString(CreateType type);
CreateType StringToCreateType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const CreateType &type);

//===--------------------------------------------------------------------===//
// Drop Types
//===--------------------------------------------------------------------===//

enum class DropType {
  INVALID = INVALID_TYPE_ID,  // invalid drop type
  DB = 1,                     // db drop type
  TABLE = 2,                  // table drop type
  INDEX = 3,                  // index drop type
  CONSTRAINT = 4,             // constraint drop type
  TRIGGER = 5,                // trigger drop type
  SCHEMA = 6,                 // trigger drop type
};
std::string DropTypeToString(DropType type);
DropType StringToDropType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const DropType &type);

template <class E>
class EnumHash {
 public:
  size_t operator()(const E &e) const {
    return std::hash<typename std::underlying_type<E>::type>()(
        static_cast<typename std::underlying_type<E>::type>(e));
  }
};

//===--------------------------------------------------------------------===//
// Language Types for UDFs
//===--------------------------------------------------------------------===//

enum class PLType {
  PL_PGSQL = 0,  // UDF language: Pl_PGSQL
  PL_C = 1       // UDF language: PL_C
};

//===--------------------------------------------------------------------===//
// Statement Types
//===--------------------------------------------------------------------===//

enum class StatementType {
  INVALID = INVALID_TYPE_ID,  // invalid statement type
  SELECT = 1,                 // select statement type
  INSERT = 3,                 // insert statement type
  UPDATE = 4,                 // update statement type
  DELETE = 5,                 // delete statement type
  CREATE = 6,                 // create statement type
  DROP = 7,                   // drop statement type
  PREPARE = 8,                // prepare statement type
  EXECUTE = 9,                // execute statement type
  RENAME = 11,                // rename statement type
  ALTER = 12,                 // alter statement type
  TRANSACTION = 13,           // transaction statement type,
  COPY = 14,                  // copy type
  ANALYZE = 15,               // analyze type
  VARIABLE_SET = 16,          // variable set statement type
  CREATE_FUNC = 17,           // create func statement type
  EXPLAIN = 18                // explain statement type
};
std::string StatementTypeToString(StatementType type);
StatementType StringToStatementType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const StatementType &type);

//===--------------------------------------------------------------------===//
// Query Types
//===--------------------------------------------------------------------===//

enum class QueryType {
  QUERY_BEGIN = 0,         // begin query
  QUERY_COMMIT = 1,        // commit query
  QUERY_ROLLBACK = 2,      // rollback query
  QUERY_CREATE_TABLE = 3,  // create query
  QUERY_CREATE_DB = 4,
  QUERY_CREATE_INDEX = 5,
  QUERY_DROP = 6,     // other queries
  QUERY_INSERT = 7,   // insert query
  QUERY_PREPARE = 8,  // prepare query
  QUERY_EXECUTE = 9,  // execute query
  QUERY_UPDATE = 10,
  QUERY_DELETE = 11,
  QUERY_RENAME = 12,
  QUERY_ALTER = 13,
  QUERY_COPY = 14,
  QUERY_ANALYZE = 15,
  QUERY_SET = 16,   // set query
  QUERY_SHOW = 17,  // show query
  QUERY_SELECT = 18,
  QUERY_OTHER = 19,
  QUERY_INVALID = 20,
  QUERY_CREATE_TRIGGER = 21,
  QUERY_CREATE_SCHEMA = 22,
  QUERY_CREATE_VIEW = 23,
  QUERY_EXPLAIN = 24
};
std::string QueryTypeToString(QueryType query_type);
QueryType StringToQueryType(std::string str);
namespace parser {
class SQLStatement;
}
QueryType StatementTypeToQueryType(StatementType stmt_type,
                                   const parser::SQLStatement *sql_stmt);
//===--------------------------------------------------------------------===//
// Scan Direction Types
//===--------------------------------------------------------------------===//

enum class ScanDirectionType {
  INVALID = INVALID_TYPE_ID,  // invalid scan direction
  FORWARD = 1,                // forward
  BACKWARD = 2                // backward
};

//===--------------------------------------------------------------------===//
// Join Types
//===--------------------------------------------------------------------===//

enum class JoinType {
  INVALID = INVALID_TYPE_ID,  // invalid join type
  LEFT = 1,                   // left
  RIGHT = 2,                  // right
  INNER = 3,                  // inner
  OUTER = 4,                  // outer
  SEMI = 5                    // IN+Subquery is SEMI
};
std::string JoinTypeToString(JoinType type);
JoinType StringToJoinType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const JoinType &type);

//===--------------------------------------------------------------------===//
// Aggregate Types
//===--------------------------------------------------------------------===//

enum class AggregateType {
  INVALID = INVALID_TYPE_ID,
  SORTED = 1,
  HASH = 2,
  PLAIN = 3  // no group-by
};
std::string AggregateTypeToString(AggregateType type);
AggregateType StringToAggregateType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const AggregateType &type);

// ------------------------------------------------------------------
// Expression Quantifier Types
// ------------------------------------------------------------------
enum class QuantifierType {
  NONE = 0,
  ANY = 1,
  ALL = 2,
};
std::string QuantifierTypeToString(QuantifierType type);
QuantifierType StringToQuantifierType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const QuantifierType &type);

//===--------------------------------------------------------------------===//
// Table Reference Types
//===--------------------------------------------------------------------===//

enum class TableReferenceType {
  INVALID = INVALID_TYPE_ID,  // invalid table reference type
  NAME = 1,                   // table name
  SELECT = 2,                 // output of select
  JOIN = 3,                   // output of join
  CROSS_PRODUCT = 4           // out of cartesian product
};
std::string TableReferenceTypeToString(TableReferenceType type);
TableReferenceType StringToTableReferenceType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const TableReferenceType &type);

//===--------------------------------------------------------------------===//
// Insert Types
//===--------------------------------------------------------------------===//

enum class InsertType {
  INVALID = INVALID_TYPE_ID,  // invalid insert type
  VALUES = 1,                 // values
  SELECT = 2                  // select
};
std::string InsertTypeToString(InsertType type);
InsertType StringToInsertType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const InsertType &type);

//===--------------------------------------------------------------------===//
// Copy Types
//===--------------------------------------------------------------------===//

enum class CopyType {
  INVALID = INVALID_TYPE_ID,  // Invalid copy type
  IMPORT_CSV,                 // Import csv data to database
  IMPORT_TSV,                 // Import tsv data to database
  EXPORT_CSV,                 // Export data to csv file
  EXPORT_STDOUT,              // Export data to std out
  EXPORT_OTHER                // Export data to other file format
};
std::string CopyTypeToString(CopyType type);
CopyType StringToCopyType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const CopyType &type);

enum class ExternalFileFormat {
  CSV,
};
std::string ExternalFileFormatToString(ExternalFileFormat format);
ExternalFileFormat StringToExternalFileFormat(const std::string &str);
std::ostream &operator<<(std::ostream &os, const ExternalFileFormat &format);

//===--------------------------------------------------------------------===//
// Payload Types
//===--------------------------------------------------------------------===//

enum class PayloadType {
  INVALID = INVALID_TYPE_ID,  // invalid message type
  CLIENT_REQUEST = 1,         // request
  CLIENT_RESPONSE = 2,        // response
  STOP = 3                    // stop loop
};
std::string PayloadTypeToString(PayloadType type);
PayloadType StringToPayloadType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const PayloadType &type);

//===--------------------------------------------------------------------===//
// Task Priority Types
//===--------------------------------------------------------------------===//

enum class TaskPriorityType {
  INVALID = INVALID_TYPE_ID,  // invalid priority
  LOW = 10,
  NORMAL = 11,
  HIGH = 12
};
std::string TaskPriorityTypeToString(TaskPriorityType type);
TaskPriorityType StringToTaskPriorityType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const TaskPriorityType &type);

//===--------------------------------------------------------------------===//
// Result Types
//===--------------------------------------------------------------------===//

enum class ResultType {
  INVALID = INVALID_TYPE_ID,  // invalid result type
  SUCCESS = 1,
  FAILURE = 2,
  ABORTED = 3,  // aborted
  NOOP = 4,     // no op
  UNKNOWN = 5,
  QUEUING = 6,
  TO_ABORT = 7,
  FAILURE_AFTER_DELETE = 8,
};
std::string ResultTypeToString(ResultType type);
ResultType StringToResultType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const ResultType &type);

//===--------------------------------------------------------------------===//
// Constraint Types
//===--------------------------------------------------------------------===//

enum class PostgresConstraintType {
  NOT_NULL, /* not standard SQL, but a lot of people * expect it */
  NOTNULL,
  DEFAULT,
  CHECK,
  PRIMARY,
  UNIQUE,
  EXCLUSION,
  FOREIGN,
  ATTR_DEFERRABLE, /* attributes for previous constraint node */
  ATTR_NOT_DEFERRABLE,
  ATTR_DEFERRED,
  ATTR_IMMEDIATE
};



ConstraintType StringToConstraintType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const ConstraintType &type);



FKConstrActionType StringToFKConstrActionType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const FKConstrActionType &type);

enum class FKConstrMatchType { SIMPLE = 0, PARTIAL = 1, FULL = 2 };

//===--------------------------------------------------------------------===//
// Set Operation Types
//===--------------------------------------------------------------------===//

enum class SetOpType {
  INVALID = INVALID_TYPE_ID,
  INTERSECT = 1,
  INTERSECT_ALL = 2,
  EXCEPT = 3,
  EXCEPT_ALL = 4
};
std::string SetOpTypeToString(SetOpType type);
SetOpType StringToSetOpType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const SetOpType &type);

// //===--------------------------------------------------------------------===//
// // Logging + Recovery Types
// //===--------------------------------------------------------------------===//

// enum class LoggingType {
//   INVALID = INVALID_TYPE_ID,
//   OFF = 1,  // turn off GC
//   ON = 2    // turn on GC
// };
// std::string LoggingTypeToString(LoggingType type);
// LoggingType StringToLoggingType(const std::string &str);
// std::ostream &operator<<(std::ostream &os, const LoggingType &type);

// enum class LogRecordType {
//   INVALID = INVALID_TYPE_ID,

//   // TransactionContext-related records
//   TRANSACTION_BEGIN = 1,
//   TRANSACTION_COMMIT = 2,

//   // Generic dml records
//   TUPLE_INSERT = 11,
//   TUPLE_DELETE = 12,
//   TUPLE_UPDATE = 13,

//   // Epoch related records
//   EPOCH_BEGIN = 21,
//   EPOCH_END = 22,
// };
// std::string LogRecordTypeToString(LogRecordType type);
// LogRecordType StringToLogRecordType(const std::string &str);
// std::ostream &operator<<(std::ostream &os, const LogRecordType &type);

// enum class CheckpointingType {
//   INVALID = INVALID_TYPE_ID,
//   OFF = 1,  // turn off GC
//   ON = 2    // turn on GC
// };
// std::string CheckpointingTypeToString(CheckpointingType type);
// CheckpointingType StringToCheckpointingType(const std::string &str);
// std::ostream &operator<<(std::ostream &os, const CheckpointingType &type);


//===--------------------------------------------------------------------===//
// Logging + Recovery Types
//===--------------------------------------------------------------------===//

// AAA_BBB
// Data stored in AAA
// Log stored in BBB
enum class LoggingType {
  INVALID = INVALID_TYPE_ID,

  // Based on write behind logging
  NVM_WBL = 1,
  SSD_WBL = 2,
  HDD_WBL = 3,

  // Based on write ahead logging
  NVM_WAL = 4,
  SSD_WAL = 5,
  HDD_WAL = 6,
};
std::string LoggingTypeToString(LoggingType type);
LoggingType StringToLoggingType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const LoggingType &type);


enum class LoggerMappingStrategyType {
  INVALID = INVALID_TYPE_ID,
  ROUND_ROBIN = 1,
  AFFINITY = 2,
  MANUAL = 3
};
std::string LoggerMappingStrategyTypeToString(LoggerMappingStrategyType type);
LoggerMappingStrategyType StringToLoggerMappingStrategyType(
    const std::string &str);
std::ostream &operator<<(std::ostream &os,
                         const LoggerMappingStrategyType &type);

enum class CheckpointType {
  INVALID = INVALID_TYPE_ID,
  NORMAL = 1,
};
std::string CheckpointTypeToString(CheckpointType type);
CheckpointType StringToCheckpointType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const CheckpointType &type);

enum ReplicationType {
  ASYNC_REPLICATION,
  SYNC_REPLICATION,
  SEMISYNC_REPLICATION
};

enum class LoggingStatusType {
  INVALID = INVALID_TYPE_ID,
  STANDBY = 1,
  RECOVERY = 2,
  LOGGING = 3,
  TERMINATE = 4,
  SLEEP = 5
};
std::string LoggingStatusTypeToString(LoggingStatusType type);
LoggingStatusType StringToLoggingStatusType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const LoggingStatusType &type);

enum class LoggerType { INVALID = INVALID_TYPE_ID, FRONTEND = 1, BACKEND = 2 };
std::string LoggerTypeToString(LoggerType type);
LoggerType StringToLoggerType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const LoggerType &type);

enum LogRecordType {
  LOGRECORD_TYPE_INVALID = INVALID_TYPE_ID,

  // Transaction-related records
  LOGRECORD_TYPE_TRANSACTION_BEGIN = 1,
  LOGRECORD_TYPE_TRANSACTION_COMMIT = 2,
  LOGRECORD_TYPE_TRANSACTION_END = 3,
  LOGRECORD_TYPE_TRANSACTION_ABORT = 4,
  LOGRECORD_TYPE_TRANSACTION_DONE = 5,

  // Generic dml records
  LOGRECORD_TYPE_TUPLE_INSERT = 11,
  LOGRECORD_TYPE_TUPLE_DELETE = 12,
  LOGRECORD_TYPE_TUPLE_UPDATE = 13,

  // DML records for Write ahead logging
  LOGRECORD_TYPE_WAL_TUPLE_INSERT = 21,
  LOGRECORD_TYPE_WAL_TUPLE_DELETE = 22,
  LOGRECORD_TYPE_WAL_TUPLE_UPDATE = 23,

  // DML records for Write behind logging
  LOGRECORD_TYPE_WBL_TUPLE_INSERT = 31,
  LOGRECORD_TYPE_WBL_TUPLE_DELETE = 32,
  LOGRECORD_TYPE_WBL_TUPLE_UPDATE = 33,

  // Record for delimiting transactions
  // includes max persistent commit_id
  LOGRECORD_TYPE_ITERATION_DELIMITER = 41,

  //DDL records 
  LOGRECORD_TYPE_DDL_TABLE_CREATE = 51,
  LOGRECORD_TYPE_DDL_TABLE_DROP = 52,
  LOGRECORD_TYPE_DDL_TABLE_ALTER = 53,
};
std::string LogRecordTypeToString(LogRecordType type);
LogRecordType StringToLogRecordType(const std::string &str);

enum class CheckpointStatus {
  INVALID = INVALID_TYPE_ID,
  STANDBY = 1,
  RECOVERY = 2,
  DONE_RECOVERY = 3,
  CHECKPOINTING = 4,
};
std::string CheckpointStatusToString(CheckpointStatus type);
CheckpointStatus StringToCheckpointStatus(const std::string &str);
std::ostream &operator<<(std::ostream &os, const CheckpointStatus &type);





std::string LayoutTypeToString(LayoutType type);
std::ostream &operator<<(std::ostream &os, const LayoutType &type);

//===--------------------------------------------------------------------===//
// Trigger Types
//===--------------------------------------------------------------------===//

//enum class TriggerType {
//  BEFORE_INSERT_ROW =
//      TRIGGER_TYPE_BEFORE | TRIGGER_TYPE_INSERT | TRIGGER_TYPE_ROW,
//  BEFORE_INSERT_STATEMENT =
//      TRIGGER_TYPE_BEFORE | TRIGGER_TYPE_INSERT | TRIGGER_TYPE_STATEMENT,
//  BEFORE_UPDATE_ROW =
//      TRIGGER_TYPE_BEFORE | TRIGGER_TYPE_UPDATE | TRIGGER_TYPE_ROW,
//  BEFORE_UPDATE_STATEMENT =
//      TRIGGER_TYPE_BEFORE | TRIGGER_TYPE_UPDATE | TRIGGER_TYPE_STATEMENT,
//  BEFORE_DELETE_ROW =
//      TRIGGER_TYPE_BEFORE | TRIGGER_TYPE_DELETE | TRIGGER_TYPE_ROW,
//  BEFORE_DELETE_STATEMENT =
//      TRIGGER_TYPE_BEFORE | TRIGGER_TYPE_DELETE | TRIGGER_TYPE_STATEMENT,
//  AFTER_INSERT_ROW =
//      TRIGGER_TYPE_AFTER | TRIGGER_TYPE_INSERT | TRIGGER_TYPE_ROW,
//  AFTER_INSERT_STATEMENT =
//      TRIGGER_TYPE_AFTER | TRIGGER_TYPE_INSERT | TRIGGER_TYPE_STATEMENT,
//  AFTER_UPDATE_ROW =
//      TRIGGER_TYPE_AFTER | TRIGGER_TYPE_UPDATE | TRIGGER_TYPE_ROW,
//  AFTER_UPDATE_STATEMENT =
//      TRIGGER_TYPE_AFTER | TRIGGER_TYPE_UPDATE | TRIGGER_TYPE_STATEMENT,
//  AFTER_DELETE_ROW =
//      TRIGGER_TYPE_AFTER | TRIGGER_TYPE_DELETE | TRIGGER_TYPE_ROW,
//  AFTER_DELETE_STATEMENT =
//      TRIGGER_TYPE_AFTER | TRIGGER_TYPE_DELETE | TRIGGER_TYPE_STATEMENT,
//  ON_COMMIT_INSERT_ROW =
//      TRIGGER_TYPE_COMMIT | TRIGGER_TYPE_INSERT | TRIGGER_TYPE_ROW,
//  ON_COMMIT_INSERT_STATEMENT =
//      TRIGGER_TYPE_COMMIT | TRIGGER_TYPE_INSERT | TRIGGER_TYPE_STATEMENT,
//  ON_COMMIT_UPDATE_ROW =
//      TRIGGER_TYPE_COMMIT | TRIGGER_TYPE_UPDATE | TRIGGER_TYPE_ROW,
//  ON_COMMIT_UPDATE_STATEMENT =
//      TRIGGER_TYPE_COMMIT | TRIGGER_TYPE_UPDATE | TRIGGER_TYPE_STATEMENT,
//  ON_COMMIT_DELETE_ROW =
//      TRIGGER_TYPE_COMMIT | TRIGGER_TYPE_DELETE | TRIGGER_TYPE_ROW,
//  ON_COMMIT_DELETE_STATEMENT =
//      TRIGGER_TYPE_COMMIT | TRIGGER_TYPE_DELETE | TRIGGER_TYPE_STATEMENT,
//};
//
//const int TRIGGER_TYPE_MAX = TRIGGER_TYPE_ROW | TRIGGER_TYPE_STATEMENT |
//                             TRIGGER_TYPE_BEFORE | TRIGGER_TYPE_AFTER |
//                             TRIGGER_TYPE_COMMIT | TRIGGER_TYPE_INSERT |
//                             TRIGGER_TYPE_DELETE | TRIGGER_TYPE_UPDATE;

//===--------------------------------------------------------------------===//
// Statistics Types
//===--------------------------------------------------------------------===//

// Statistics Collection Type
// Disable or enable
// TODO: This should probably be a collection level and not a boolean
// (enable/disable)
enum class StatsType {
  // Disable statistics collection
  INVALID = INVALID_TYPE_ID,
  // Enable statistics collection
  ENABLE = 1,
};

enum class MetricType {
  // Metric type is invalid
  INVALID = INVALID_TYPE_ID,
  // Metric to count a number
  COUNTER = 1,
  // Access information, e.g., # tuples read, inserted, updated, deleted
  ACCESS = 2,
  // Life time of a object
  LIFETIME = 3,
  // Statistics for a specific database
  DATABASE = 4,
  // Statistics for a specific table
  TABLE = 5,
  // Statistics for a specific index
  INDEX = 6,
  // Latency of transactions
  LATENCY = 7,
  // Timestamp, e.g., creation time of a table/index
  TEMPORAL = 8,
  // Statistics for a specific table
  QUERY = 9,
  // Statistics for CPU
  PROCESSOR = 10,
};

// All builtin operators we currently support
enum class OperatorId : uint32_t {
  Negation = 0,
  Abs,
  Add,
  Sub,
  Mul,
  Div,
  Mod,
  LogicalAnd,
  LogicalOr,
  Ascii,
  Chr,
  Concat,
  Substr,
  CharLength,
  OctetLength,
  Length,
  Repeat,
  Replace,
  LTrim,
  RTrim,
  BTrim,
  Trim,
  Sqrt,
  Ceil,
  Round,
  DatePart,
  Floor,
  DateTrunc,
  Like,
  Now,

  // Add more operators here, before the last "Invalid" entry
  Invalid
};
std::string OperatorIdToString(OperatorId op_id);

enum class OnError : uint32_t { ReturnNull, Exception };

static const int INVALID_FILE_DESCRIPTOR = -1;

// ------------------------------------------------------------------
// Tuple serialization formats
// ------------------------------------------------------------------

enum class TupleSerializationFormat { NATIVE = 0, DR = 1 };

// ------------------------------------------------------------------
// Entity types
// ------------------------------------------------------------------

enum class EntityType {
  INVALID = INVALID_TYPE_ID,
  TABLE = 1,
  SCHEMA = 2,
  INDEX = 3,
  VIEW = 4,
  PREPARED_STATEMENT = 5,
};
std::string EntityTypeToString(EntityType type);
EntityType StringToEntityType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const EntityType &type);

//===--------------------------------------------------------------------===//
// Type definitions.
//===--------------------------------------------------------------------===//

typedef size_t hash_t;

typedef uint32_t uint;





static const uint MAX_OID = std::numeric_limits<uint>::max() - 1;

#define NULL_OID MAX_OID









static const cid_t INVALID_EID = 0;

static const cid_t MAX_EID = std::numeric_limits<eid_t>::max();

// For epoch
//====================wjh
static const size_t EPOCH_LENGTH = 40;
static const size_t ACTIVE_COUNT = 40;
static const size_t MT_PS_COUNT = 40;
//static const size_t ACTIVE_COUNT = 1000;
//static const size_t MT_PS_COUNT = 1000;
static const size_t MEM_POOL_COUNT = 100;
//static const bool FOR_NVM = false;
//static const bool FOR_MEM = true;
static const bool FOR_NVM = true;
static const bool FOR_MEM = false;
static const bool INDEX_LOG_OPEN = false;

static const size_t chpt_COUNT = 1000000000;//十亿条 2000仓

static const IndexType INDEXTYPE = IndexType::MASSTREE;
static const uint INDEXPARM = 0;

//static const IndexType indexType = IndexType::BWTREE;
//static const uint indexParm = 1;

static const size_t TASK_NO = 1000;
static const size_t THD_NO = 1;
static const size_t CHECKPOINT_TIME = 3000;
//====================wjh

//====================wjh

static const size_t per_tile_group_count = 1000;
static const size_t RETRYTIME = 5;
static const size_t Every_RETRY_TIME = 5;
// For threads
extern size_t CONNECTION_THREAD_COUNT;
extern size_t LOGGING_THREAD_COUNT;
extern size_t GC_THREAD_COUNT;
extern size_t EPOCH_THREAD_COUNT;


/** @define Align size to 8 bytes. */
#define ALIGN8(len) (((len) + 8 - 1) & ~((size_t)(8 - 1)))
#define GetBytes4(x) (((uintptr_t)(x)) & 0xffffffff)
#define GetBytes8(x) ((uintptr_t)(x))
//===--------------------------------------------------------------------===//
// TupleMetadata
//===--------------------------------------------------------------------===//
struct TupleMetadata {
  uint table_id = 0;
  uint tile_group_id = 0;
  uint tuple_slot_id = 0;
  cid_t tuple_end_cid = 0;
};

//===--------------------------------------------------------------------===//
// Column Bitmap
//===--------------------------------------------------------------------===//
static const size_t max_col_count = 128;
typedef std::bitset<max_col_count> ColBitmap;

//===--------------------------------------------------------------------===//
// read-write set
//===--------------------------------------------------------------------===//

// this enum is to identify the operations performed by the transaction.
enum class RWType {
  INVALID,
  READ,
  READ_OWN,  // select for update
  UPDATE,
  INSERT,
  DELETE,
  INS_DEL,  // delete after insert.
};
std::string RWTypeToString(RWType type);
RWType StringToRWType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const RWType &type);

// ItemPointer -> type
typedef tbb::concurrent_unordered_map<ItemPointer, RWType, ItemPointerHasher,
                                      ItemPointerComparator>
    ReadWriteSet;

//typedef std::unordered_map<ItemPointer, RWType, ItemPointerHasher,
//    ItemPointerComparator>
//    ReadWriteSet;

//typedef std::unordered_set<ItemPointer, ItemPointerHasher,
//                                      ItemPointerComparator>
//    WriteSet;
typedef tbb::concurrent_unordered_set<ItemPointer, ItemPointerHasher,
                                      ItemPointerComparator>
    WriteSet;

// this enum is to identify why the version should be GC'd.
enum class GCVersionType {
  INVALID,
  COMMIT_UPDATE,   // a version that is updated during txn commit.
  COMMIT_DELETE,   // a version that is deleted during txn commit.
  COMMIT_INS_DEL,  // a version that is inserted and deleted during txn commit.
  ABORT_UPDATE,    // a version that is updated during txn abort.
  ABORT_DELETE,    // a version that is deleted during txn abort.
  ABORT_INSERT,    // a version that is inserted during txn abort.
  ABORT_INS_DEL,   // a version that is inserted and deleted during txn commit.
};
std::string GCVersionTypeToString(GCVersionType type);
GCVersionType StringToGCVersionType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const GCVersionType &type);

// block -> offset -> type
typedef std::unordered_map<uint, std::unordered_map<uint, GCVersionType>>
    GCSet;

enum class DDLType {
  INVALID,
  CREATE,
  DROP,
};
typedef tbb::concurrent_vector<std::tuple<uint, uint, uint, DDLType>>
    CreateDropSet;

//typedef std::vector<std::tuple<uint, uint, uint, DDLType>>
//    CreateDropSet;
typedef std::vector<std::tuple<uint, uint, uint>> GCObjectSet;

//===--------------------------------------------------------------------===//
// File Handle
//===--------------------------------------------------------------------===//
struct FileHandle {
  // FILE pointer
  FILE *file = nullptr;

  // File descriptor
  int fd;

  // Size of the file
  size_t size;

  FileHandle() : file(nullptr), fd(INVALID_FILE_DESCRIPTOR), size(0) {}

  FileHandle(FILE *file, int fd, size_t size)
      : file(file), fd(fd), size(size) {}
};
extern FileHandle INVALID_FILE_HANDLE;

//===--------------------------------------------------------------------===//
// Transformers
//===--------------------------------------------------------------------===//



std::string TypeIdArrayToString(const std::vector<TypeId> &types);
std::vector<TypeId> StringToTypeArray(const std::string &types);


std::string SqlStateErrorCodeToString(SqlStateErrorCode code);

//namespace expression {
//class AbstractExpression;
//}

/**
 * @brief Generic specification of a projection target:
 *        < DEST_column_id , expression >
 */

//namespace planner {
//struct DerivedAttribute;
//}
//
//typedef std::pair<uint, const DerivedAttribute> Target;

//typedef std::vector<Target> TargetList;

/**
 * @brief Generic specification of a direct map:
 *        < NEW_col_id , <tuple_index (left or right tuple), OLD_col_id>    >
 */
typedef std::pair<uint, std::pair<uint, uint>> DirectMap;

typedef std::vector<DirectMap> DirectMapList;

//===--------------------------------------------------------------------===//
// Optimizer
//===--------------------------------------------------------------------===//

enum class PropertyType {
  INVALID = INVALID_TYPE_ID,
  COLUMNS,
  DISTINCT,
  SORT,
  LIMIT,
};
std::string PropertyTypeToString(PropertyType type);
PropertyType StringToPropertyType(const std::string &str);
std::ostream &operator<<(std::ostream &os, const PropertyType &type);

enum class RuleType : uint32_t {
  // Transformation rules (logical -> logical)
  INNER_JOIN_COMMUTE = 0,
  INNER_JOIN_ASSOCIATE,

  // Don't move this one
  LogicalPhysicalDelimiter,

  // Implementation rules (logical -> physical)
  GET_TO_DUMMY_SCAN,
  GET_TO_SEQ_SCAN,
  GET_TO_INDEX_SCAN,
  QUERY_DERIVED_GET_TO_PHYSICAL,
  EXTERNAL_FILE_GET_TO_PHYSICAL,
  DELETE_TO_PHYSICAL,
  UPDATE_TO_PHYSICAL,
  INSERT_TO_PHYSICAL,
  INSERT_SELECT_TO_PHYSICAL,
  AGGREGATE_TO_HASH_AGGREGATE,
  AGGREGATE_TO_PLAIN_AGGREGATE,
  INNER_JOIN_TO_NL_JOIN,
  INNER_JOIN_TO_HASH_JOIN,
  IMPLEMENT_DISTINCT,
  IMPLEMENT_LIMIT,
  EXPORT_EXTERNAL_FILE_TO_PHYSICAL,

  // Don't move this one
  RewriteDelimiter,

  // Rewrite rules (logical -> logical)
  PUSH_FILTER_THROUGH_JOIN,
  COMBINE_CONSECUTIVE_FILTER,
  EMBED_FILTER_INTO_GET,
  MARK_JOIN_GET_TO_INNER_JOIN,
  MARK_JOIN_INNER_JOIN_TO_INNER_JOIN,
  MARK_JOIN_FILTER_TO_INNER_JOIN,
  PULL_FILTER_THROUGH_MARK_JOIN,
  PULL_FILTER_THROUGH_AGGREGATION,

  // Place holder to generate number of rules compile time
  NUM_RULES

};

//class AbstractExpression;
//class ExprHasher;
//class ExprEqualCmp;

// Augment abstract expression with a table alias set
//struct AnnotatedExpression {
//  AnnotatedExpression(std::shared_ptr<AbstractExpression> i_expr,
//                      std::unordered_set<std::string> &i_set)
//      : expr(i_expr), table_alias_set(i_set) {}
//  AnnotatedExpression(const AnnotatedExpression &mt_expr)
//      : expr(mt_expr.expr), table_alias_set(mt_expr.table_alias_set) {}
//  std::shared_ptr<AbstractExpression> expr;
//  std::unordered_set<std::string> table_alias_set;
//};

//typedef std::vector<AnnotatedExpression> MultiTablePredicates;
//typedef std::unordered_map<
//    std::string, std::vector<std::shared_ptr<AbstractExpression>>>
//    SingleTablePredicatesMap;

// TODO(boweic): use raw ptr
// Mapping of Expression -> Column Offset created by operator
//typedef std::unordered_map<AbstractExpression *, unsigned,
//                           ExprHasher, ExprEqualCmp>
//    ExprMap;
//// Used in optimizer to speed up expression comparsion
//typedef std::unordered_set<AbstractExpression *,
//                           ExprHasher, ExprEqualCmp>
//    ExprSet;



//===--------------------------------------------------------------------===//
// Wire protocol typedefs
//===--------------------------------------------------------------------===//
#define SOCKET_BUFFER_SIZE 8192

/* byte type */
typedef unsigned char uchar;

/* type for buffer of bytes */
typedef std::vector<uchar> ByteBuf;

//===--------------------------------------------------------------------===//
// Packet Manager: ProcessResult
//===--------------------------------------------------------------------===//
enum class ProcessResult {
  COMPLETE,
  TERMINATE,
  PROCESSING,
  MORE_DATA_REQUIRED,
  NEED_SSL_HANDSHAKE,
};

enum class NetworkProtocolType {
  POSTGRES_JDBC,
  POSTGRES_PSQL,
};

enum class SSLLevel {
  SSL_DISABLE = 0,
  SSL_PREFER = 1,
  SSL_VERIIFY = 2,
};

// Eigen/Matrix types used in brain
// TODO(saatvik): Generalize Eigen utilities across all types
typedef std::vector<float> vector_t;
typedef std::vector<std::vector<float>> matrix_t;
//typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
//    matrix_eig;
//typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> vector_eig;
