-- Main database tables
CREATE TABLE IF NOT EXISTS `libraries` (
  `library_id` INTEGER PRIMARY KEY,
  `code` TEXT NOT NULL,
  `name` TEXT NOT NULL,
  `source` TEXT NOT NULL,
	`total_reads` INT NOT NULL
);

CREATE TABLE IF NOT EXISTS `sequences` (
  `sid` INTEGER PRIMARY KEY,
  `seq` TEXT UNIQUE NOT NULL
);

CREATE TABLE IF NOT EXISTS `reads` (
  `sid` INT NOT NULL,
  `library_id` INT NOT NULL,
  `reads` INT NOT NULL
);
 
CREATE UNIQUE INDEX IF NOT EXISTS `sid` ON `reads` (`sid`,`library_id`);
CREATE INDEX IF NOT EXISTS `library_id` ON `reads` (`library_id`);

-- Features tables
CREATE TABLE IF NOT EXISTS `feature_categories` (
  `category_id` INT PRIMARY KEY,
  `category` TEXT NOT NULL DEFAULT '' UNIQUE,
  `source` TEXT NOT NULL DEFAULT ''
);

INSERT INTO `feature_categories` (`category_id`, `category`, `source`) VALUES(1, 'gene', '');
INSERT INTO `feature_categories` (`category_id`, `category`, `source`) VALUES(2, 'mRNA', '');
INSERT INTO `feature_categories` (`category_id`, `category`, `source`) VALUES(3, 'five_prime_UTR', '');
INSERT INTO `feature_categories` (`category_id`, `category`, `source`) VALUES(4, 'exon', '');
INSERT INTO `feature_categories` (`category_id`, `category`, `source`) VALUES(5, 'three_prime_UTR', '');
INSERT INTO `feature_categories` (`category_id`, `category`, `source`) VALUES(6, 'CDS', '');
INSERT INTO `feature_categories` (`category_id`, `category`, `source`) VALUES(7, 'DNA_transposon', '');
INSERT INTO `feature_categories` (`category_id`, `category`, `source`) VALUES(8, 'LTR_retrotransposon', '');
INSERT INTO `feature_categories` (`category_id`, `category`, `source`) VALUES(9, 'non-LTR_retrotransposon', '');
INSERT INTO `feature_categories` (`category_id`, `category`, `source`) VALUES(10, 'other_repeat', '');
INSERT INTO `feature_categories` (`category_id`, `category`, `source`) VALUES(11, 'unknown_repeat', '');

CREATE TABLE IF NOT EXISTS `features` (
  `feature_id` INT PRIMARY KEY,
  `category_id` INT NOT NULL,
  `ref_seq` TEXT NOT NULL,
  `start` INT NOT NULL,
  `end` INT NOT NULL,
  `strand` INT NOT NULL,
  `accession` TEXT NOT NULL,
  `sub_category` TEXT NOT NULL,
  `description` TEXT NOT NULL,
  `parent_id` INT NOT NULL
);

CREATE INDEX IF NOT EXISTS `accession` ON `features` (`accession`);
CREATE INDEX IF NOT EXISTS `start` ON `features` (`start`);
CREATE INDEX IF NOT EXISTS `end` ON `features` (`end`);
CREATE INDEX IF NOT EXISTS `ref_seq` ON `features` (`ref_seq`);
