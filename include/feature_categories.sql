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
