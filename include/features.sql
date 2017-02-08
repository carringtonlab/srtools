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
